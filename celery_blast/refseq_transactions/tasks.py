from os import getcwd, chdir, mkdir, remove
from os.path import isdir, isfile
from subprocess import Popen, PIPE as subPIPE, STDOUT as subSTDOUT
import pandas as pd
import math

from blast_project.py_django_db_services import get_database_by_id, update_blast_database_with_task_result_model
from .py_services import write_blastdatabase_snakemake_configfile,\
    get_ftp_paths_and_taxids_from_summary_file, write_progress_database_transactions,\
    get_current_progress_database_transactions, get_bdb_summary_table_name

from time import sleep
from celery import shared_task
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder
#logger for celery worker instances
logger = get_task_logger(__name__)

@shared_task()
def write_alias_file(alias_filename,available_databases):
    try:
        #logger.info('starting to write database alias file : {}'.format(alias_filename))
        alias_file = open(alias_filename,'w')
        database_title = '.'.join(alias_filename.split("/")[3].split(".")[0:2])
        alias_file.write("TITLE {}\n".format(database_title))
        alias_file.write("DBLIST")
        for database_file in available_databases:
            alias_file.write(" \""+database_file.split("/")[3]+"\"")
        alias_file.write("\n")
        alias_file.close()
        return 0
    except Exception as e:
        logger.warning("error during creation of aliasfile with exception : {}".format(e))
        raise Exception("error during creation of aliasfile with exception : {}".format(e))

@shared_task()
def format_blast_databases(path_to_database,chunks, progress_recorder):
    try:
        errorlist=[]
        available_databases=[]

        progress = 75
        # to 75 %
        steps = 100 / len(chunks) / 4
        logger.info("starting to format database chunks")
        progress_recorder.set_progress(progress, 100, "starting to format database chunks")

        for chunk in chunks:
            database = path_to_database + "database_chunk_{}.faa".format(chunk)
            taxonomic_mapfile = path_to_database + "acc_taxmap_file_{}.table".format(chunk)
            proc = Popen(["makeblastdb","-in",database,
                                     '-dbtype','prot',
                                     '-taxid_map',taxonomic_mapfile,
                                     '-parse_seqids',
                                     '-out',database,'-blastdb_version','5'])
            returncode = proc.wait(timeout=6000)
            if(returncode != 0):
                logger.warning("ERROR during database creation of chunk {}".format(chunk))
                errorlist.append(chunk)
            else:
                available_databases.append(database)
                logger.info("\tprocessed database chunk {}".format(chunk))
                progress += steps
                progress_recorder.set_progress(progress, 100, "formatted chunk {}".format(chunk))

        progress_recorder.set_progress(99, 100, "formatted chunk {}".format(chunk))
        return available_databases,errorlist

    except Exception as e:
        logger.warning('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))
        raise Exception('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))

@shared_task()
def create_chunks_of_databases(df, path_to_database, progress_recorder):
    try:
        total_formatted = 0
        chunk = 1
        chunks = []

        progress = 50
        #to 75 %
        steps = 100 / len(df['ftp_path']) / 4
        logger.info("starting to concatenate available fasta assemblies")
        progress_recorder.set_progress(progress, 100, "concatenate available fasta assemblies")

        iteration_steps = 500
        if round(len(df['ftp_path']) / iteration_steps) > 35:
            iteration_steps = math.ceil( len(df['ftp_path']) / 35)

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
                    remove(path_to_database+assembly_name)

                    for line in lines:
                        if line[0] == ">":
                            split = line.split(" ")
                            header = ' '.join(split[1:])
                            acc_id = split[0].split(">")[1] + "_" + fasta_header
                            #accession header max length for makeblastdb is 50
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
    except Exception as e:
        logger.warning('couldnt create database chunks with exception : {}'.format(e))
        raise Exception('couldnt create database chunks with exception : {}'.format(e))
    return chunks


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
    
'''
@shared_task()
def download_wget_ftp_paths(path_to_database,dictionary_ftp_paths_taxids,progress_recorder):
    try:
        error_log = open(path_to_database+'download_error.log', 'w')
        transform_ftp_path = lambda file: file.split('/')[-1].rstrip(file[-3:])

        #progress_steps = round(100 / (len(dictionary_ftp_paths_taxids.keys()) * 2), 3)

        downloaded_files = {}
        #database_id for format_available_databases

        progress_steps = 100/((len(dictionary_ftp_paths_taxids.keys()))*2)
        progress = 0
        progress_recorder.set_progress(progress, 100, "starting wget processes")

        for file in dictionary_ftp_paths_taxids.keys():

            gunzip_output = transform_ftp_path(file)
            if(isfile(path_to_database + gunzip_output) == True):
                downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                logger.info('file {} exists, scipping download'.format(gunzip_output))
                progress += progress_steps
                progress_recorder.set_progress(progress, 100, "downloaded {}".format(gunzip_output))
            elif(isfile(path_to_database + gunzip_output) == False):
                #0 ... 9 attempts
                for attempt in range(10):
                    try:


                        proc = Popen('wget -qO- {} | gzip -d > {}'.format(file, path_to_database+gunzip_output), shell=True)
                        returncode = proc.wait(timeout=300)  # 66 Minutes
                        if(returncode != 0):
                            raise Exception
                        #downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                        #logger.info('downloaded : {} returncode : {}'.format(path_to_database + gunzip_output,returncode))

                    #catch exception raised if the subproccess failed e.g. gzip failure due to invalid download
                    except Exception as e:
                        logger.warning("download exception : {}".format(e))
                        error_log.write("{} {}\n".format(file, attempt))
                        logger.warning('next download attempt of file : {} with attempt : {}'.format(file,attempt))
                        if(attempt == 9):
                            error_log.write('couldnt download: {} '.format(file))
                            logger.warning('couldnt download: {} after 10 attempts'.format(file))
                            if(isfile(path_to_database + gunzip_output) == True):
                                remove(path_to_database + gunzip_output)
                                logger.warning("removed empty file: {}".format(gunzip_output))

                            progress += progress_steps
                            progress_recorder.set_progress(progress, 100, "failed trying to download {}".format(file))
                            #break
                    else:
                        logger.info('downloaded : {} returncode : {}'.format(path_to_database + gunzip_output,returncode))
                        downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                        progress += progress_steps
                        progress_recorder.set_progress(progress, 100, "downloaded {}".format(gunzip_output))
                        break

        error_log.close()
        return downloaded_files
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
        db_path = media/databases/2/
        progress_recorder = CeleryProgressBar Class
        dict = get_ftp_paths_and_taxids_from_summary_file(database_id)
        dict = download_blast_databases_based_on_summary_file(db_path,dict,progress_recorder)
        df = pd.DictToDf(dict)
        chunks = create_chunks_of_databases(db_path,df,progress_recorder)
        available_db_chunks, errorlist = format_blast_databases(db_path,chunks,progress_recorder)
        
'''
#TODO documentation task 1
@shared_task(bind=True)
def download_blast_databases_based_on_summary_file(self, database_id):
    try:
        progress_recorder = ProgressRecorder(self)
        logger.info('starting downloading task')
        working_dir = getcwd()
        logger.info('working dir : {}'.format(working_dir))
        path_to_database = 'media/databases/' + str(database_id) + '/'

        progress_recorder.set_progress(0,100,"started downloading")

        #loads the desired database summary table into an dictionary with the 'ftp_path' column as key
        #and 'taxid' as value
        dictionary_ftp_paths_taxids = get_ftp_paths_and_taxids_from_summary_file(database_id)

        logger.info('trying to update blast database with current task')
        update_blast_database_with_task_result_model(database_id, str(self.request.id))

        #download
        dictionary_ftp_paths_taxids = download_wget_ftp_paths(path_to_database,dictionary_ftp_paths_taxids,progress_recorder)

        #build pandas dataframe from dictionary that is returned by the download function
        pandas_ftp_paths_taxids_df = pd.DataFrame(dictionary_ftp_paths_taxids.items(), columns=['ftp_path', 'taxid'])
        available_database_chunks = create_chunks_of_databases(pandas_ftp_paths_taxids_df,path_to_database,progress_recorder)

        #execute makeblastdb
        databases = format_blast_databases(path_to_database,available_database_chunks,progress_recorder)

        #database_files = format_available_databases(path_to_database,dictionary_ftp_paths_taxids,progress_recorder)
        progress_recorder.set_progress(99, 100, "writing alias file")
        database_pandas_table_name = get_bdb_summary_table_name(database_id)
        alias_filename = 'media/databases/' + str(database_id) + '/' + database_pandas_table_name+'.database.pal'

        logger.info('starting to write database alias file : {}'.format(alias_filename))

        returncode = write_alias_file(alias_filename,databases[0])
        '''
        alias_file = open(alias_filename,'w')
        alias_file.write("TITLE {}\n".format(database_pandas_table_name+'.database'))
        alias_file.write("DBLIST")
        for database_file in database_files:
            alias_file.write(" \""+database_file+"\"")
        alias_file.write("\n")
        alias_file.close()
        '''

        progress_recorder.set_progress(100, 100, "writing alias file")
        return returncode
    except Exception as e:
        logger.warning('couldnt perform task because of exception : {}'.format(e))
        raise Exception('couldnt perform compute ftp path task with exception {}'.format(e))


#TODO documentation
#task_track_started --> no need for celery-progress?
@shared_task(bind=True)
def download_blast_databases(self, database_id):


    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100,"started process")

    try:
        update_blast_database_with_task_result_model(database_id, str(self.request.id))
    except Exception as e:
        logger.warning("couldnt save taskresult into blastdatabase ... ")
        raise Exception("coulndt save taskresult into blastdatabase ...")

    #check if blast summary file is available
    blast_database_summary_table = get_database_by_id(database_id).get_pandas_table_name()
    logger.info('load blastdatabase: {} from database'.format(blast_database_summary_table))
    bdb_summary_table_path = 'media/databases/' + str(database_id) + '/' + blast_database_summary_table
    if (isfile(bdb_summary_table_path) == False):
        logger.warning('there is no database table with path : {}', format(bdb_summary_table_path))
        raise Exception(
            "couldnt download assembly files, there is no summary table with path : {}".format(bdb_summary_table_path))

    # summary file and ProgressRecord
    try:
        write_blastdatabase_snakemake_configfile(database_id, str(self.request.id))
    except Exception as e:
        logger.warning("couldnt write snakemake configuration file")
        raise Exception("couldnt write snakemake configuration file")


    progress_recorder.set_progress(20, 100, "starting snakemake process")

    try:
        snakemake_working_dir = 'media/databases/' + str(database_id) + '/'
        snakemake_config_file = 'media/databases/' + str(database_id) + '/snakefile_config'
        snakefile_dir = 'static/snakefiles/assembly_download/Snakefile'
        #cmd to copy
        #snakemake --snakefile 'static/snakefiles/assembly_download/Snakefile' --cores 2 --configfile 'media/databases/12/snakefile_config' --directory 'media/databases/12/' --wms-monitor 'http://172.23.0.5:5000' --latency-wait 10 --dry-run
        snakemake_process = Popen(['snakemake',
                                   '--snakefile',snakefile_dir,
                                   '--wms-monitor','http://172.23.0.6:5000',
                                   '--cores','1',
                                   '--configfile',snakemake_config_file,
                                   '--directory',snakemake_working_dir,
                                   '--latency-wait','10',
                                   '--quiet',
                                   '--keep-incomplete'], shell=False, stdout=subPIPE, stderr=subSTDOUT)
        logger.info('waiting for popen snakemake instance {} to finish'.format(snakemake_process.pid))
        progress_recorder.set_progress(40, 100, "waiting for snakemake to finish")
        returncode = snakemake_process.wait(timeout=4000) # 66 Minutes

        if (returncode != 0):
            logger.warning('subprocess Popen snakemake instance resulted in an error!')
            raise Exception('Popen hasnt succeeded ...')

        logger.info('returncode : {}'.format(returncode))
        logger.info('blast database download and formatting completed')
        progress_recorder.set_progress(100, 100, "complete")
        return returncode
    except Exception as e:
        raise Exception("couldnt download assembly files, error during execution of snakemake, with exception : ".format(e))

''' download_refseq_assembly_summary_file
    
    this function downloads the current assembly_summary_refseq.txt file into 
    the media directory /databases/refseq_summary_file/ if this directory does not exists the function
    creates the directory
    
    it uses subprocess.Popen to invoke wget, stdout (e.g. percentage is printed out inside celery worker console)
    
    :returns download_file
        :type str
'''
@shared_task()
def download_refseq_assembly_summary_file():
    try:
        refseq_url = "ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
        current_working_directory = getcwd() #/blast/reciprocal_blast
        path_to_assembly_file_location = current_working_directory + '/media/databases/refseq_summary_file/'
        timeout=300

        logger.info('setup filepath parameter:\n\t cwd : {} \n\t path_to_assembly_file_location : {}'
                    .format(current_working_directory,path_to_assembly_file_location))

        if(isdir(path_to_assembly_file_location) == False):
            logger.warning('path_to_assembly_file_location : {} does not exists, trying to create it with mkdir ...')
            mkdir(path_to_assembly_file_location)

        path_to_assembly_file = path_to_assembly_file_location+'assembly_summary_refseq.txt'
        if(isfile(path_to_assembly_file)):
            logger.warning('assembly_summary_refseq.txt exists deleting it in order to download a newer version')
            remove(path_to_assembly_file)

        # invoke wget program

        wget_process = Popen(['wget', refseq_url,'-q','-O',path_to_assembly_file], shell=False,stdout=subPIPE, stderr=subSTDOUT)
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

        return refseq_url.split('/')[-1]
    except Exception as e:
        raise Exception("couldn't download assembly_summary_refseq.txt file : {}".format(e))


'''
#TODO documentation
@shared_task()
def format_available_databases(path_to_database, database_files_taxids_dict, progress_recorder):
    try:
        progress_steps = 100 / ((len(database_files_taxids_dict.keys())) * 2)
        progress = 50
        progress_recorder.set_progress(progress, 100, "starting makeblastdb processes")
        for assembly_file in database_files_taxids_dict.keys():
            database = path_to_database + assembly_file
            taxid = database_files_taxids_dict[assembly_file]
            cmd = "makeblastdb -in {} -dbtype prot -out {} -taxid {}".format(database,database,taxid)

            if (isfile(database+'.pdb') == True):
                logger.info('file {} exists, scipping formatting'.format(database))
            else:
                proc = Popen("makeblastdb -in {} -dbtype prot -out {} -taxid {} -parse_seqids"
                             .format(database,database,taxid), shell=True)
                returncode = proc.wait(timeout=600)#10*60 Sec.
                logger.info('processed database with makeblastdb : {}'.format(cmd))
                if(returncode != 0):
                    logger.warning("something went wrong during creation of database: {}".format(database))

            progress += progress_steps
            progress_recorder.set_progress(progress, 100, "starting makeblastdb processes")

        progress_recorder.set_progress(99, 100, "ended database formatting")
        return database_files_taxids_dict.keys()
    except Exception as e:
        logger.warning('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))
        raise Exception('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))
'''
