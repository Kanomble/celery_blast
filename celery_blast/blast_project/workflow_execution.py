import uuid
from dataclasses import dataclass
from typing import Optional

from django.db import transaction
from django_celery_results.models import TaskResult

from celery_blast.resource_governance import enforce_user_workflow_quota

from .models import BlastProject, RemoteBlastProject, WorkflowLifecycle


ACTIVE_TASK_STATUSES = {'PENDING', 'STARTED', 'PROGRESS', 'RETRY'}
STARTABLE_STATES = {
    WorkflowLifecycle.STARTABLE,
    WorkflowLifecycle.SUCCESSFUL,
    WorkflowLifecycle.FAILED,
    WorkflowLifecycle.CANCELLED,
}

TASK_STATUS_TO_LIFECYCLE = {
    'PENDING': WorkflowLifecycle.QUEUED,
    'STARTED': WorkflowLifecycle.RUNNING,
    'PROGRESS': WorkflowLifecycle.RUNNING,
    'RETRY': WorkflowLifecycle.RETRYING,
    'SUCCESS': WorkflowLifecycle.SUCCESSFUL,
    'FAILURE': WorkflowLifecycle.FAILED,
    'REVOKED': WorkflowLifecycle.CANCELLED,
}


@dataclass(frozen=True)
class WorkflowStartResult:
    started: bool
    project_id: int
    task_id: Optional[str]
    lifecycle_state: str


def lifecycle_from_task_status(status):
    if not status:
        return WorkflowLifecycle.STARTABLE
    return TASK_STATUS_TO_LIFECYCLE.get(status, WorkflowLifecycle.RUNNING)


def is_active_task_status(status):
    return status in ACTIVE_TASK_STATUSES


def _project_config(remote=False):
    if remote:
        return {
            'model': RemoteBlastProject,
            'task_field': 'r_project_execution_snakemake_task',
            'state_field': 'r_project_workflow_state',
            'token_field': 'r_project_workflow_execution_token',
            'detail_route': 'remote_project_details',
            'project_type': 'remote',
            'user_field': 'r_project_user',
        }
    return {
        'model': BlastProject,
        'task_field': 'project_execution_snakemake_task',
        'state_field': 'project_workflow_state',
        'token_field': 'project_workflow_execution_token',
        'detail_route': 'project_details',
        'project_type': 'local',
        'user_field': 'project_user',
    }


def _refresh_lifecycle_from_task(project, config):
    task_result = getattr(project, config['task_field'])
    state_field = config['state_field']
    token_field = config['token_field']
    current_state = getattr(project, state_field)

    if task_result is None:
        if current_state not in STARTABLE_STATES:
            setattr(project, state_field, WorkflowLifecycle.STARTABLE)
            setattr(project, token_field, None)
        return getattr(project, state_field), None

    lifecycle_state = lifecycle_from_task_status(task_result.status)
    if lifecycle_state != current_state:
        setattr(project, state_field, lifecycle_state)
    if getattr(project, token_field) is None and task_result.task_id:
        setattr(project, token_field, task_result.task_id)
    return lifecycle_state, task_result


def _enqueue_default_task(project_id, execution_token, remote=False):
    if remote:
        from .tasks import execute_remote_reciprocal_blast_project

        execute_remote_reciprocal_blast_project.apply_async(
            args=(project_id, execution_token),
            task_id=execution_token,
        )
    else:
        from .tasks import execute_reciprocal_blast_project

        execute_reciprocal_blast_project.apply_async(
            args=(project_id, execution_token),
            task_id=execution_token,
        )


def start_reciprocal_workflow(project_id, remote=False, task_sender=None):
    config = _project_config(remote=remote)
    task_sender = task_sender or (lambda pk, token: _enqueue_default_task(pk, token, remote=remote))

    with transaction.atomic():
        project = (
            config['model'].objects
            .select_for_update()
            .get(id=project_id)
        )
        lifecycle_state, task_result = _refresh_lifecycle_from_task(project, config)
        if (
            lifecycle_state not in STARTABLE_STATES
            or (task_result is not None and is_active_task_status(task_result.status))
        ):
            project.save(update_fields=[config['state_field'], config['token_field']])
            return WorkflowStartResult(
                started=False,
                project_id=project.id,
                task_id=getattr(task_result, 'task_id', None),
                lifecycle_state=lifecycle_state,
            )

        enforce_user_workflow_quota(getattr(project, config['user_field']))

        execution_token = uuid.uuid4().hex
        task_result = TaskResult.objects.create(task_id=execution_token, status='PENDING')
        setattr(project, config['task_field'], task_result)
        setattr(project, config['state_field'], WorkflowLifecycle.QUEUED)
        setattr(project, config['token_field'], execution_token)
        project.save(update_fields=[config['task_field'], config['state_field'], config['token_field']])

        transaction.on_commit(lambda: task_sender(project.id, execution_token))

        return WorkflowStartResult(
            started=True,
            project_id=project.id,
            task_id=execution_token,
            lifecycle_state=WorkflowLifecycle.QUEUED,
        )


def claim_workflow_execution(project_id, execution_token, project_type):
    config = _project_config(remote=(project_type == 'remote'))
    with transaction.atomic():
        project = config['model'].objects.select_for_update().get(id=project_id)
        if getattr(project, config['token_field']) != execution_token:
            return False
        if getattr(project, config['state_field']) != WorkflowLifecycle.QUEUED:
            return False
        task_result = TaskResult.objects.get(task_id=execution_token)
        setattr(project, config['task_field'], task_result)
        setattr(project, config['state_field'], WorkflowLifecycle.RUNNING)
        project.save(update_fields=[config['task_field'], config['state_field']])
        return True


def finish_workflow_execution(project_id, execution_token, project_type, successful):
    config = _project_config(remote=(project_type == 'remote'))
    with transaction.atomic():
        project = config['model'].objects.select_for_update().get(id=project_id)
        if getattr(project, config['token_field']) != execution_token:
            return False
        setattr(
            project,
            config['state_field'],
            WorkflowLifecycle.SUCCESSFUL if successful else WorkflowLifecycle.FAILED,
        )
        project.save(update_fields=[config['state_field']])
        return True
