from subprocess import Popen, SubprocessError, TimeoutExpired

from celery import shared_task
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder
from django.conf import settings

from .py_django_db_services import update_one_way_blast_project_with_task_result_model, \
    update_one_way_remote_blast_project_with_task_result_model

# logger for celery worker instances
logger = get_task_logger(__name__)


@shared_task(bind=True)
def execute_one_way_blast_project(self, project_id):
    logger.info('snakemake one way blast workflow execution')
    logger.info('setting panoptes ip to: {}'.format(settings.PANOPTES_IP))
    snakemake_working_dir = settings.ONE_WAY_BLAST_PROJECT_DIR + str(project_id) + '/'
    snakemake_config_file = settings.ONE_WAY_BLAST_PROJECT_DIR + str(project_id) + '/snakefile_config'
    snakefile_dir = 'static/snakefiles/one_way_blast/Snakefile'

    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100, "started process")

    try:
        update_one_way_blast_project_with_task_result_model(project_id, str(self.request.id))
    except Exception as e:
        logger.warning('couldnt update blastproject with exception : {}'.format(e))
        raise Exception('couldnt update blastproject with exception : {}'.format(e))
    try:
        logger.info('trying to start snakemake one way BLAST workflow')
        # snakemake --snakefile 'static/snakefiles/reciprocal_blast/Snakefile' --cores 1 --configfile 'media/blast_projects/4/snakefile_config' --directory 'media/blast_projects/4/'
        # '--wms-monitor', settings.PANOPTES_IP,

        one_way_blast_snakemake = Popen(
            ['snakemake',
             '--snakefile', snakefile_dir,
             '--cores', '1',
             '--configfile', snakemake_config_file,
             '--directory', snakemake_working_dir,
             '--keep-incomplete'], shell=False)
        progress_recorder.set_progress(50, 100, "executed snakemake")

        logger.info('waiting for popen instance {} to finish with timeout set to {}'.format(one_way_blast_snakemake.pid,
                                                                                            settings.SUBPROCESS_TIME_LIMIT))
        returncode = one_way_blast_snakemake.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)
        if returncode > 0:
            logger.warning('received a negative returncode : {} ... '.format(returncode))
            raise Exception
        logger.info('returncode : {}'.format(returncode))
        progress_recorder.set_progress(100, 100, "SUCCESS")
        return returncode

    except TimeoutExpired as e:
        logger.info('timeout expired ... trying to kill process {}'.format(one_way_blast_snakemake.pid))
        one_way_blast_snakemake.kill()

        raise Exception(
            'exception during waiting for popen instance : {} \n\t returncode of popen.wait : {}'.format(e, returncode))
    except SubprocessError as e:
        logger.info('subprocess throwed exception: {}'.format(e))
        raise Exception(
            'exception occured during invokation of:\n\t reciprocal blast snakefile : {}'.format(e))


@shared_task(bind=True)
def execute_one_way_remote_blast_project(self, project_id):
    logger.info('snakemake one way remote blast workflow execution')
    logger.info('setting panoptes ip to: {}'.format(settings.PANOPTES_IP))
    snakemake_working_dir = settings.ONE_WAY_BLAST_PROJECT_DIR + '/remote_searches/' + str(project_id) + '/'
    snakemake_config_file = settings.ONE_WAY_BLAST_PROJECT_DIR + '/remote_searches/' + str(
        project_id) + '/snakefile_config'
    snakefile_dir = 'static/snakefiles/one_way_blast/remote_searches/Snakefile'

    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100, "started process")

    try:
        update_one_way_remote_blast_project_with_task_result_model(project_id, str(self.request.id))
    except Exception as e:
        logger.warning('couldnt update blastproject with exception : {}'.format(e))
        raise Exception('couldnt update blastproject with exception : {}'.format(e))
    try:
        logger.info('trying to start snakemake one way remote BLAST workflow')
        # snakemake --snakefile 'static/snakefiles/one_way_blast/remote_searches/Snakefile' --cores 1 --configfile
        # 'media/one_way_blast/remote_searches/6/snakefile_config' --directory 'media/one_way_blast/remote_searches/6/'
        #             '--wms-monitor', settings.PANOPTES_IP,
        one_way_remote_blast_snakemake = Popen(
            ['snakemake',
             '--snakefile', snakefile_dir,
             '--cores', '1',
             '--configfile', snakemake_config_file,
             '--directory', snakemake_working_dir,
             '--keep-incomplete'], shell=False)
        progress_recorder.set_progress(50, 100, "executed snakemake")

        logger.info(
            'waiting for popen instance {} to finish with timeout set to {}'.format(one_way_remote_blast_snakemake.pid,
                                                                                    settings.SUBPROCESS_TIME_LIMIT))
        returncode = one_way_remote_blast_snakemake.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)
        if returncode > 0:
            logger.warning('received a negative returncode {} ... '.format(returncode))
            raise Exception
        logger.info('returncode : {}'.format(returncode))
        progress_recorder.set_progress(100, 100, "SUCCESS")
        return returncode

    except TimeoutExpired as e:
        logger.info('timeout expired ... trying to kill process {}'.format(one_way_remote_blast_snakemake.pid))
        one_way_remote_blast_snakemake.kill()
        raise Exception(
            'exception during waiting for popen instance : {} \n\t returncode of popen.wait > 0'.format(e))

    except SubprocessError as e:
        logger.info('subprocess throwed exception: {}'.format(e))
        raise Exception(
            'exception occured during invokation of:\n\t reciprocal blast snakefile : {}'.format(e))
