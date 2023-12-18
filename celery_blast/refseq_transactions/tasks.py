import math
import os
from os import getcwd, mkdir, remove, chmod
from os.path import isdir, isfile
from subprocess import Popen, SubprocessError, TimeoutExpired

import pandas as pd
from blast_project.py_django_db_services import update_blast_database_with_task_result_model, update_assembly_entries_in_database
from celery import shared_task
from celery.exceptions import SoftTimeLimitExceeded
from celery.utils.log import get_task_logger
from celery_blast.settings import BLAST_DATABASE_DIR, REFSEQ_URL, REFSEQ_ASSEMBLY_FILE, GENBANK_ASSEMBLY_FILE, GENBANK_URL
from celery_progress.backend import ProgressRecorder
from django.conf import settings
from .py_services import get_ftp_paths_and_taxids_from_summary_file, get_bdb_summary_table_name, update_blast_database_table

# logger for celery worker instances
logger = get_task_logger(__name__)

''' download_refseq_assembly_summary_file

    This function downloads the current assembly_summary_refseq.txt file or the assembly_summary_genbank.txt file into 
    the media directory /databases/refseq_summary_file/ if this directory does not exists the function
    creates the directory.
    WARNING: This function uses the default BLAST_DATABASE_DIR filepath, overwriting BLAST_DATABASE_DIR does not
    affect the function.
    
    it uses subprocess.Popen to invoke wget, stdout (e.g. percentage is printed out inside celery worker console)

    :returns returncode of Popen
        :type str
'''


@shared_task()
def download_refseq_assembly_summary(summary_file:str):
    try:
        current_working_directory = getcwd()  # /blast/reciprocal_blast
        if summary_file == "RefSeq":
            path_to_assembly_file_location = REFSEQ_ASSEMBLY_FILE
            refseq_url = REFSEQ_URL

        elif summary_file == "GenBank":
            path_to_assembly_file_location = GENBANK_ASSEMBLY_FILE
            refseq_url = GENBANK_URL
        else:
            logger.warning("wrong summary_file selected: {}".format(summary_file))
            logger.warning("specify RefSeq or GenBank in the relevant form.html document ...")
            raise Exception("wrong summary_file selected: {}".format(summary_file))
        timeout = 300

        logger.info('setup filepath parameter:\n\t cwd : {} \n\t path_to_assembly_file_location : {}'
                    .format(current_working_directory, path_to_assembly_file_location))

        if (isdir(path_to_assembly_file_location) == False):
            logger.warning('path_to_assembly_file_location : {} does not exists, trying to create it with mkdir ...')
            mkdir(path_to_assembly_file_location)

        if summary_file == "RefSeq":
            path_to_assembly_file = REFSEQ_ASSEMBLY_FILE + 'assembly_summary_refseq.txt'
        elif summary_file == "GenBank":
            path_to_assembly_file = GENBANK_ASSEMBLY_FILE + 'assembly_summary_genbank.txt'

        if (isfile(path_to_assembly_file)):
            logger.warning('assembly_summary_refseq.txt/assembly_summary_genbank.txt exists deleting it in order to download a newer version')
            remove(path_to_assembly_file)

        # invoke wget program
        logger.info('creating popen process')
        logger.info("refseq_url: {}, path_to_assembly_file: {}".format(refseq_url, path_to_assembly_file))
        wget_process = Popen(['wget', refseq_url, '-q', '-O', path_to_assembly_file], shell=False)
        # communicate with subprocess : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.communicate
        # wait for process to terminate and set returncode attribute
        logger.info(
            'waiting for popen instance {} to finish with timeout set to {}'
                .format(wget_process.pid, timeout))
        returncode = wget_process.wait(timeout=timeout)

        if (returncode != 0):
            logger.warning('subprocess Popen refseq assembly file download process resulted in an error!')
            raise Exception('Popen hasnt succeeded ...')

        logger.info('returncode : {}'.format(returncode))

        logger.info('download completed')

        return returncode

    except SoftTimeLimitExceeded:
        if 'path_to_assembly_file' in locals():
            if isfile(path_to_assembly_file):
                os.remove(path_to_assembly_file)
        raise Exception("ERROR couldn't download assembly_summary_refseq.txt file due to soft time limit")

    except TimeoutExpired as e:
        wget_process.kill()
        if 'path_to_assembly_file' in locals():
            if isfile(path_to_assembly_file):
                os.remove(path_to_assembly_file)

        raise Exception(
            "ERROR couldn't download assembly_summary_refseq.txt file due Popen call time limit : {}".format(e))

    except SubprocessError as e:
        wget_process.kill()
        if 'path_to_assembly_file' in locals():
            if isfile(path_to_assembly_file):
                os.remove(path_to_assembly_file)
        raise Exception(
            "ERROR couldn't download assembly_summary_refseq.txt file due Popen call exception : {}".format(e))

    except Exception as e:
        raise Exception("couldn't download assembly_summary_refseq.txt file : {}".format(e))


f'''write_alias_file

    This function produces a .pal file, which lists all fasta files belonging to the dedicated database.
    
    :param alias_filename - Name of the .pal file.
        :type str
        
    :param available_databases - List of all database_chunks_integer.fasta
        :type list
        
    :returns 0 on success
        :type int
'''


@shared_task()
def write_alias_file(alias_filename: str, available_databases: list) -> int:
    try:
        logger.info('starting to write database alias file : {}'.format(alias_filename))
        with open(alias_filename, 'w') as alias_file:
            database_title = '.'.join(alias_filename.split("/")[3].split(".")[0:2])
            alias_file.write("TITLE {}\n".format(database_title))
            alias_file.write("DBLIST")
            for database_file in available_databases:
                alias_file.write(" \"" + database_file.split("/")[3] + "\"")
            alias_file.write("\n")
            alias_file.close()
        return 0
    except SoftTimeLimitExceeded:
        if 'alias_file' in locals():
            if isfile(alias_file):
                os.remove(alias_file)
        raise Exception("ERROR writing alias file due to soft time limit")

    except Exception as e:
        logger.warning("error during creation of aliasfile with exception : {}".format(e))
        if 'alias_file' in locals():
            if isfile(alias_file):
                os.remove(alias_file)
        raise Exception("error during creation of aliasfile with exception : {}".format(e))


'''format_blast_databases

    This function executes the makeblastdb command in an subprocess Popen instance.
    Makeblastdb is used to format fasta files into BLAST databases. 
    The function operates on all database chunks that are provided via a list of integers. 
    E.g. [0,1,2] will produce three database files: database_chunk_1.faa ...
    
    :param path_to_database - Filepath, e.g. media/databases/1
        :type str
    
    :param chunks - List of integers.
        :type list
    
    :param progress_recorder - Used to update the celery task.
        :type ProgressRecorder class
    
    :returns available databases - List of successfully formatted database names
        :type list
    
    :returns errorlist - List of databases with errors during the formatting procedure.
        :type list
'''


@shared_task()
def format_blast_databases(path_to_database: str, chunks: list, progress_recorder: ProgressRecorder) -> tuple:
    try:
        errorlist = []
        available_databases = []

        progress = 75
        # to 75 %
        steps = 100 / len(chunks) / 4
        logger.info("starting to format database chunks")
        progress_recorder.set_progress(progress, 100, "starting to format database chunks")

        for chunk in chunks:
            database = path_to_database + "database_chunk_{}.faa".format(chunk)
            taxonomic_mapfile = path_to_database + "acc_taxmap_file_{}.table".format(chunk)
            proc = Popen(["makeblastdb", "-in", database,
                          '-dbtype', 'prot',
                          '-taxid_map', taxonomic_mapfile,
                          '-parse_seqids',
                          '-out', database, '-blastdb_version', '5'])
            returncode = proc.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)
            if (returncode != 0):
                logger.warning("ERROR during database creation of chunk {}".format(chunk))
                errorlist.append(chunk)
            else:
                available_databases.append(database)
                logger.info("\tprocessed database chunk {}".format(chunk))
                progress += steps
                progress_recorder.set_progress(progress, 100, "formatted chunk {}".format(chunk))

        progress_recorder.set_progress(99, 100, "formatted chunk {}".format(chunk))
        return available_databases, errorlist

    except TimeoutExpired as e:
        logger.warning("ERROR in makeblastdb for database chunk: {} - process timed out".format(database))
        if 'proc' in locals():
            proc.kill()
        raise Exception(
            "ERROR in makeblastdb for database chunk: {} with exception: {} - process timed out".format(database, e))

    except SubprocessError as e:
        logger.warning("ERROR in makeblastdb for database chunk: {}".format(database))
        if 'proc' in locals():
            proc.kill()
        raise Exception("ERROR in makeblastdb for database chunk: {} with exception: {}".format(database, e))

    except SoftTimeLimitExceeded:
        raise Exception("ERROR in formatting blast databases with makeblastdb soft time limit exceeded ...")
    except Exception as e:
        logger.warning('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))
        raise Exception('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))


'''create_chunks_of_databases
    This function combines each 500 fasta files to one database chunk that can gets formatted to 
    BLAST databases by makeblastdb. The provided pandas dataframe (dict[str,int]) is based on the output
    of download_wget_ftp_paths.
    
    :param df - based on dict[str, int] with str=filepath and int=taxid
        :type pd.DataFrame
    
    :param path_to_database
        :type str
    
    :param progress_recorder
        :type ProgressRecorder
        
    :returns chunks - list of integers that define the database chunks
        :type list - list[int]
'''


@shared_task()
def create_chunks_of_databases(df: pd.DataFrame, path_to_database: str, progress_recorder: ProgressRecorder) -> list:
    try:
        chunk_logfile = open(path_to_database+'create_chunks_of_databases.log','w')

        total_formatted = 0
        chunk = 1
        chunks = []

        progress = 50
        # to 75 %
        steps = 100 / len(df['ftp_path']) / 4
        logger.info("starting to concatenate available fasta assemblies")
        progress_recorder.set_progress(progress, 100, "concatenate available fasta assemblies")

        iteration_steps = 500  # chunks --> 500 fasta files build "one" database chunk
        # the following two lines ensure that there are no more than 35 different databases in the pal file,
        # thus the iteration_steps will be adjusted to fit the 35 database files ...
        if round(len(df['ftp_path']) / iteration_steps) > 35:
            iteration_steps = math.ceil(len(df['ftp_path']) / 35)

        logger.info("set iteration steps to : {}".format(iteration_steps))
        while total_formatted < len(df['ftp_path']):
            iteration_end = total_formatted + iteration_steps
            if (iteration_end > len(df['ftp_path'])):
                iteration_end = len(df['ftp_path'])

            logger.info("looping from {} to {}".format(total_formatted, iteration_end))

            if (isfile(path_to_database + "acc_taxmap_file_{}.table".format(chunk)) == True and
                    isfile(path_to_database + "database_chunk_{}.faa".format(chunk)) == True):
                logger.info("Skipping chunk creation, database exists")
                # previous format process of assemblies : df['ftp_path'][total_formatted:iteration_end]
                total_formatted = iteration_end
                for it in range(iteration_steps):
                    progress += steps
                progress_recorder.set_progress(progress, 100, "skipped available chunk")
                logger.info("total formatted assembly_files: {}".format(total_formatted))
                chunks.append(chunk)
                chunk += 1
            else:
                acc_to_tax = open(path_to_database + "acc_taxmap_file_{}.table".format(chunk), 'w')
                database_chunk = open(path_to_database + "database_chunk_{}.faa".format(chunk), 'w')

                for ftp_path, taxid in zip(df['ftp_path'][total_formatted:iteration_end],
                                           df['taxid'][total_formatted:iteration_end]):
                    assembly_name = ftp_path
                    logger.info('\tworking with {}'.format(assembly_name))

                    fasta_header = '_'.join(assembly_name.split('_')[0:3])
                    genome_file = open(path_to_database + assembly_name, 'r')

                    lines = genome_file.readlines()
                    genome_file.close()

                    try:
                        remove(path_to_database + assembly_name)
                    except:
                        chunk_logfile.write("failed to remove file: {}\n".format(path_to_database + assembly_name))
                        logger.warning("[-] couldnt remove file: {}".format(path_to_database + assembly_name))

                    for line in lines:
                        # transformes accession id and adds additional information
                        if line[0] == ">":
                            split = line.split(" ")
                            header = ' '.join(split[1:])
                            acc_id = split[0].split(">")[1] + "_" + fasta_header
                            # accession header max length for makeblastdb is 50
                            acc_to_tax.write(acc_id[0:49] + "\t" + str(taxid) + "\n")
                            line = '>' + acc_id[0:49] + ' ' + header
                        database_chunk.write(line)

                    total_formatted += 1
                    progress += steps
                    progress_recorder.set_progress(progress, 100, "processed chunk {}".format(chunk))

                logger.info("\tdone writing database chunk {}".format(chunk))
                database_chunk.close()
                chunks.append(chunk)
                chunk += 1
                acc_to_tax.close()

        chunk_logfile.close()
        return chunks
    except SoftTimeLimitExceeded:
        logger.warning("ERROR soft time limit exceeded for creation of database chunks")
        raise Exception("ERROR soft time limit exceeded for creation of database chunks")
    except Exception as e:
        logger.warning('couldnt create database chunks with exception : {}'.format(e))
        raise Exception('couldnt create database chunks with exception : {}'.format(e))


''' download_wget_ftp_paths
    
    Downloads assembly files that reside in the Database summary table by their corresponding ftp_path. 
    Assemblies are directly decompressed with gzip. This function gets executed by download_blast_databases_based_on_summary_file.
    The returned dictionary can be used to concatenate the fasta files into database chunks.
    
    :param path_to_database
        type: str - e.g. media/databases/2/
    :param dictionary_ftp_paths_taxids
        type: dict - dictionary with ftp_paths as keys and taxids as values
    :param progress_recorder
        type: class ProgressRecorder - celery progress bars class
        
    :returns downloaded_files
        type: dict - dictionary of downloaded filenames and their corresponding taxids as values
    :returns errorlist
        type: list[str]
    
'''


@shared_task()
def download_wget_ftp_paths(path_to_database: str, dictionary_ftp_paths_taxids: dict,
                            progress_recorder: ProgressRecorder) -> dict:
    try:
        error_log = open(path_to_database + 'download_error.log', 'w')
        transform_ftp_path = lambda file: file.split('/')[-1].rstrip(file[-3:])

        # progress_steps = round(100 / (len(dictionary_ftp_paths_taxids.keys()) * 2), 3)

        downloaded_files = {}
        # database_id for format_available_databases

        progress_steps = 100 / ((len(dictionary_ftp_paths_taxids.keys())) * 2)
        progress = 0
        progress_recorder.set_progress(progress, 100, "starting wget processes")
        errorlist = []
        for file in dictionary_ftp_paths_taxids.keys():

            gunzip_output = transform_ftp_path(file)
            if (isfile(path_to_database + gunzip_output) == True):
                downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                logger.info('file {} exists, scipping download'.format(gunzip_output))
                progress += progress_steps
                progress_recorder.set_progress(progress, 100, "downloaded {}".format(gunzip_output))
            elif (isfile(path_to_database + gunzip_output) == False):
                # 0 ... 9 attempts
                for attempt in range(10):
                    try:

                        proc = Popen('wget -qO- {} | gzip -d > {}'.format(file, path_to_database + gunzip_output),
                                     shell=True)
                        returncode = proc.wait(timeout=300)  # 66 Minutes
                        if (returncode != 0):
                            raise Exception
                        # downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                        # logger.info('downloaded : {} returncode : {}'.format(path_to_database + gunzip_output,returncode))
                        chmod(path_to_database + gunzip_output, 0o777)
                    # catch exception raised if the subproccess failed e.g. gzip failure due to invalid download
                    except Exception as e:
                        logger.warning("download exception : {}".format(e))
                        error_log.write("{} {}\n".format(file, attempt))
                        logger.warning('next download attempt of file : {} with attempt : {}'.format(file, attempt))
                        if (attempt == 9):
                            error_log.write('couldnt download: {} '.format(file))
                            logger.warning('couldnt download: {} after 10 attempts'.format(file))
                            if (isfile(path_to_database + gunzip_output) == True):
                                remove(path_to_database + gunzip_output)  # removes partially downloaded files
                                logger.warning("removed empty file: {}".format(gunzip_output))

                            progress += progress_steps
                            progress_recorder.set_progress(progress, 100, "failed trying to download {}".format(file))
                            errorlist.append(file)
                            # break
                    # fills the returned dictionary with successfully downloaded genome files
                    else:
                        logger.info(
                            'downloaded : {} returncode : {}'.format(path_to_database + gunzip_output, returncode))
                        downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                        progress += progress_steps
                        progress_recorder.set_progress(progress, 100, "downloaded {}".format(gunzip_output))
                        break

        error_log.close()
        return downloaded_files, errorlist
    except Exception as e:
        logger.warning('couldnt download assemblies with exception : {}'.format(e))
        raise Exception('couldnt download assemblies with exception : {}'.format(e))


''' download_blast_databases_based_on_summary_file
    
    Main download and database formatting function, that manages temporary output.
    
    :param self
        :type asynchronous task - used for the TaskResult model
    :param database_id
        :type int -  
    
    summary of important functions and variables inside this task:
        db_path = BLAST_DATABASE_DIR + 2/
        progress_recorder = CeleryProgressBar Class
        dict = get_ftp_paths_and_taxids_from_summary_file(database_id)
        dict = download_blast_databases_based_on_summary_file(db_path,dict,progress_recorder)
        df = pd.DictToDf(dict)
        chunks = create_chunks_of_databases(db_path,df,progress_recorder)
        available_db_chunks, errorlist = format_blast_databases(db_path,chunks,progress_recorder)
        
'''


@shared_task(bind=True)
def download_blast_databases_based_on_summary_file(self, database_id):
    try:
        progress_recorder = ProgressRecorder(self)
        logger.info('starting downloading task')
        working_dir = getcwd()
        logger.info('working dir : {}'.format(working_dir))
        path_to_database = BLAST_DATABASE_DIR + str(database_id) + '/'

        progress_recorder.set_progress(0, 100, "started downloading")

        # loads the desired database summary table into an dictionary with the 'ftp_path' column as key
        # and 'taxid' as value
        dictionary_ftp_paths_taxids = get_ftp_paths_and_taxids_from_summary_file(database_id)

        logger.info('trying to update blast database with current task')
        update_blast_database_with_task_result_model(database_id, str(self.request.id))

        # download
        dictionary_ftp_paths_taxids, errorlist = download_wget_ftp_paths(path_to_database, dictionary_ftp_paths_taxids,
                                                              progress_recorder)
        logger.info('delete unused assemblies from blast database table')
        database_pandas_table_name = get_bdb_summary_table_name(database_id)
        path_to_pandas_table = path_to_database + database_pandas_table_name
        download_error_logfile = path_to_database + "update_blast_database_table.log"
        returncode = update_blast_database_table(errorlist, path_to_pandas_table, download_error_logfile)



        # build pandas dataframe from dictionary that is returned by the download function
        pandas_ftp_paths_taxids_df = pd.DataFrame(dictionary_ftp_paths_taxids.items(), columns=['ftp_path', 'taxid'])
        available_database_chunks = create_chunks_of_databases(pandas_ftp_paths_taxids_df, path_to_database,
                                                               progress_recorder)

        # execute makeblastdb
        databases = format_blast_databases(path_to_database, available_database_chunks, progress_recorder)

        # database_files = format_available_databases(path_to_database,dictionary_ftp_paths_taxids,progress_recorder)
        progress_recorder.set_progress(99, 100, "writing alias file")
        alias_filename = BLAST_DATABASE_DIR + str(database_id) + '/' + database_pandas_table_name + '.database.pal'

        logger.info('starting to write database alias file : {}'.format(alias_filename))

        returncode = write_alias_file(alias_filename, databases[0])



        returncode = update_assembly_entries_in_database(database_id)
        progress_recorder.set_progress(100, 100, "writing alias file")
        return returncode
    except Exception as e:
        logger.warning('couldnt perform task with exception : {}'.format(e))
        raise Exception('couldnt perform task with exception {}'.format(e))
