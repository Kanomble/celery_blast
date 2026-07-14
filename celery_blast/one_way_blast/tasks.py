import os

from celery import shared_task
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder
from django.conf import settings

from celery_blast.processes import ExternalCommandError, ExternalCommandTimeout, run_external_command

from .py_django_db_services import update_one_way_blast_project_with_task_result_model, \
    update_one_way_remote_blast_project_with_task_result_model

# logger for celery worker instances
logger = get_task_logger(__name__)


def resolve_project_path(path):
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(settings.BASE_DIR, path))


def build_snakemake_environment():
    env = os.environ.copy()
    pythonpath = env.get('PYTHONPATH')
    if pythonpath:
        env['PYTHONPATH'] = settings.BASE_DIR + os.pathsep + pythonpath
    else:
        env['PYTHONPATH'] = settings.BASE_DIR
    return env


def execute_snakemake_workflow(project_id, working_dir, config_file, snakefile_dir, task_result_updater, task_id,
                               progress_recorder):
    try:
        task_result_updater(project_id, task_id)
    except Exception as e:
        logger.warning('couldnt update blastproject with exception : {}'.format(e))
        raise Exception('couldnt update blastproject with exception : {}'.format(e))

    logger.info('trying to start snakemake workflow')
    command = [
        'snakemake',
        '--snakefile', resolve_project_path(snakefile_dir),
        '--cores', '1',
        '--configfile', resolve_project_path(config_file),
        '--directory', resolve_project_path(working_dir),
        '--keep-incomplete',
    ]
    result = run_external_command(
        command,
        timeout=settings.SUBPROCESS_TIME_LIMIT,
        shell=False,
        logger=logger,
        check=True,
        env=build_snakemake_environment(),
    )
    progress_recorder.set_progress(100, 100, "SUCCESS")
    return result.returncode


@shared_task(bind=True)
def execute_one_way_blast_project(self, project_id):
    logger.info('snakemake one way blast workflow execution')
    snakemake_working_dir = settings.ONE_WAY_BLAST_PROJECT_DIR + str(project_id) + '/'
    snakemake_config_file = settings.ONE_WAY_BLAST_PROJECT_DIR + str(project_id) + '/snakefile_config'
    snakefile_dir = 'static/snakefiles/one_way_blast/Snakefile'

    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100, "started process")

    try:
        progress_recorder.set_progress(50, 100, "executed snakemake")
        return execute_snakemake_workflow(
            project_id,
            snakemake_working_dir,
            snakemake_config_file,
            snakefile_dir,
            update_one_way_blast_project_with_task_result_model,
            str(self.request.id),
            progress_recorder,
        )
    except ExternalCommandTimeout as e:
        raise Exception('exception during waiting for popen instance : {}'.format(e))
    except ExternalCommandError as e:
        raise Exception(
            'exception occured during invokation of:\n\t one way blast snakefile : {}'.format(e))


@shared_task(bind=True)
def execute_one_way_remote_blast_project(self, project_id):
    logger.info('snakemake one way remote blast workflow execution')
    snakemake_working_dir = settings.ONE_WAY_BLAST_PROJECT_DIR + '/remote_searches/' + str(project_id) + '/'
    snakemake_config_file = settings.ONE_WAY_BLAST_PROJECT_DIR + '/remote_searches/' + str(
        project_id) + '/snakefile_config'
    snakefile_dir = 'static/snakefiles/one_way_blast/remote_searches/Snakefile'

    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100, "started process")

    try:
        progress_recorder.set_progress(50, 100, "executed snakemake")
        return execute_snakemake_workflow(
            project_id,
            snakemake_working_dir,
            snakemake_config_file,
            snakefile_dir,
            update_one_way_remote_blast_project_with_task_result_model,
            str(self.request.id),
            progress_recorder,
        )
    except ExternalCommandTimeout as e:
        raise Exception('exception during waiting for popen instance : {}'.format(e))
    except ExternalCommandError as e:
        raise Exception(
            'exception occured during invokation of:\n\t one way remote blast snakefile : {}'.format(e))
