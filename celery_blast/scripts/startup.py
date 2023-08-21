'''SCRIPT THAT GETS EXECUTED DURING SERVER STARTUP
1. Adding assembly level models to database
2. Downloading and extracting taxonomy database
3. Downloading and extracting CDD database
'''

from refseq_transactions.models import AssemblyLevels
from external_tools.models import DomainDatabase
from blast_project.py_services import check_domain_database_status

from os.path import isfile, isdir
from os import remove, getcwd, mkdir, listdir
import subprocess
import psutil
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, TAXDB_URL, CDD_DATABASE_URL

def run():
    if len(AssemblyLevels.objects.all()) != 4:
        print("INFO:Inserting Assembly Levels")
        AssemblyLevels.objects.all().delete()
        contig = AssemblyLevels(assembly_level='Contig')
        chromosome = AssemblyLevels(assembly_level='Chromosome')
        complete = AssemblyLevels(assembly_level='Complete Genome')
        scaffold = AssemblyLevels(assembly_level='Scaffold')
        contig.save()
        chromosome.save()
        complete.save()
        scaffold.save()

    if len(DomainDatabase.objects.all()) != 1:
        DomainDatabase.objects.all().delete()
        domain_database_model = DomainDatabase(domain_database_loaded=False)
        domain_database_model.save()
    else:
        domain_database = DomainDatabase.objects.all()[0]
        if domain_database.domain_database_loaded == False:
            check_domain_database_status()

    if isfile(BLAST_DATABASE_DIR + 'taxdb.btd') and isfile(BLAST_DATABASE_DIR + 'taxdb.bti'):
        print("INFO:TAXONOMY DATABASE IS LOADED")

    else:
        print("INFO:NO TAXONOMY DATABASE")
        if isfile(BLAST_DATABASE_DIR + "taxdb.tar.gz"):
            remove(BLAST_DATABASE_DIR + "taxdb.tar.gz")
        if isfile(BLAST_DATABASE_DIR + "taxdb.tar"):
            remove(BLAST_DATABASE_DIR + "taxdb.tar")
        print("INFO:STARTING TO DOWNLOAD TAXONOMY DATABASE")

        try:
            taxdb_ftp_path = TAXDB_URL
            current_working_directory = getcwd()  # /blast/reciprocal_blast
            path_to_taxdb_location = current_working_directory + BLAST_DATABASE_DIR
            path_to_taxdb_location = path_to_taxdb_location + 'taxdb.tar.gz'

            proc = subprocess.Popen(["wget", taxdb_ftp_path, "-q", "-O", path_to_taxdb_location], shell=False)
            returncode = proc.wait(timeout=600)
            if returncode != 0:
                raise subprocess.SubprocessError
            print("INFO:EXTRACTING TAXONOMY DB")
            database_dir = "/blast/reciprocal_blast/"+BLAST_DATABASE_DIR
            proc = subprocess.Popen(["tar", "-zxvf", path_to_taxdb_location, "-C", database_dir], shell=False)
            returncode = proc.wait(timeout=600)
            if returncode != 0:
                raise subprocess.SubprocessError
            print("INFO:DONE DOWNLOADING TAXONOMY DATABASE")

        except subprocess.TimeoutExpired as e:
            print("ERROR:TIMEOUT EXPIRED DURING DOWNLOAD OF TAXONOMY DATABASE: {}".format(e))
            print(
                "INFO:IF YOU HAVE NO STABLE INTERNET CONNECTION TRY TO RESTART THE WEBSERVER ONCE YOU HAVE A BETTER CONNECTION")
            print("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")
            if 'proc' in locals():
                pid = proc.pid
                parent = psutil.Process(pid)

                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()
            else:
                print("WARNING: CHECK FOR UNFINISHED PROCESSES OR RESTART THE WEB-SERVER")
        except subprocess.SubprocessError as e:
            print("ERROR:WGET RESULTED IN AN ERROR: {}".format(e))
            print("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")

            if 'proc' in locals():
                pid = proc.pid
                parent = psutil.Process(pid)

                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()
            else:
                print("WARNING: CHECK FOR UNFINISHED PROCESSES OR RESTART THE WEB-SERVER")