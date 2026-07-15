from django.conf import settings
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _


ACTIVE_TASK_STATUSES = {'PENDING', 'STARTED', 'PROGRESS', 'RETRY'}


class WorkflowQuotaExceeded(ValidationError):
    pass


def positive_int_setting(name, default=1):
    value = int(getattr(settings, name, default))
    return max(1, value)


def configured_limit(name, default):
    return positive_int_setting(name, default)


def workflow_job_cores():
    return positive_int_setting('CATHI_WORKFLOW_JOB_CORES', 1)


def effective_blast_threads():
    return max(
        1,
        min(
            workflow_job_cores(),
            positive_int_setting('CATHI_MAX_BLAST_THREADS', 1),
        ),
    )


def capped_positive_int(value, setting_name, default):
    value = max(1, int(value))
    return min(value, positive_int_setting(setting_name, default))


def capped_blast_threads(value):
    return min(max(1, int(value)), effective_blast_threads())


def capped_blast_settings_values(num_threads, num_alignments, max_target_seqs=None, max_hsps=None):
    values = {
        'num_threads': capped_blast_threads(num_threads),
        'num_alignments': capped_positive_int(num_alignments, 'CATHI_MAX_NUM_ALIGNMENTS', 10000),
    }
    if max_target_seqs is not None:
        values['max_target_seqs'] = capped_positive_int(max_target_seqs, 'CATHI_MAX_TARGET_SEQS', 10000)
    if max_hsps is not None:
        values['max_hsps'] = capped_positive_int(max_hsps, 'CATHI_MAX_HSPS', 500)
    return values


def validate_uploaded_file_size(uploaded_file):
    max_bytes = positive_int_setting('CATHI_MAX_UPLOAD_BYTES', 100 * 1024 * 1024)
    size = getattr(uploaded_file, 'size', 0) or 0
    if size > max_bytes:
        raise ValidationError(
            _('Uploaded file exceeds the configured %(limit)s byte limit.'),
            code='file_too_large',
            params={'limit': max_bytes},
        )
    return uploaded_file


def resource_budget_log_line():
    return (
        f"Resource budget: snakemake_cores={workflow_job_cores()}, "
        f"blast_threads={effective_blast_threads()}, "
        f"max_upload_bytes={positive_int_setting('CATHI_MAX_UPLOAD_BYTES', 100 * 1024 * 1024)}"
    )


def _reciprocal_active_count(user):
    from blast_project.models import BlastProject, RemoteBlastProject, WorkflowLifecycle

    active_states = [
        WorkflowLifecycle.QUEUED,
        WorkflowLifecycle.RUNNING,
        WorkflowLifecycle.RETRYING,
    ]
    return (
        BlastProject.objects.filter(project_user=user, project_workflow_state__in=active_states).count()
        + RemoteBlastProject.objects.filter(r_project_user=user, r_project_workflow_state__in=active_states).count()
    )


def _one_way_active_count(user):
    from one_way_blast.models import OneWayBlastProject, OneWayRemoteBlastProject

    return (
        OneWayBlastProject.objects.filter(
            project_user=user,
            project_execution_task_result__status__in=ACTIVE_TASK_STATUSES,
        ).count()
        + OneWayRemoteBlastProject.objects.filter(
            r_project_user=user,
            r_project_execution_task_result__status__in=ACTIVE_TASK_STATUSES,
        ).count()
    )


def active_workflow_count_for_user(user):
    if not user or not getattr(user, 'is_authenticated', False):
        return 0
    return _reciprocal_active_count(user) + _one_way_active_count(user)


def enforce_user_workflow_quota(user):
    max_active = positive_int_setting('CATHI_MAX_ACTIVE_WORKFLOWS_PER_USER', 1)
    if active_workflow_count_for_user(user) >= max_active:
        raise WorkflowQuotaExceeded(
            _('You already have the maximum number of active workflows.'),
            code='workflow_quota_exceeded',
        )
