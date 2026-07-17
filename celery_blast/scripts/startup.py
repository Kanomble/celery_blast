'''SCRIPT THAT GETS EXECUTED DURING SERVER STARTUP
1. Adding assembly level models to database
2. Downloading and extracting taxonomy database
3. Checks if CDD database folder exists.
'''

from refseq_transactions.models import AssemblyLevels
from external_tools.models import DomainDatabase
from blast_project.py_services import check_domain_database_status

from os.path import isfile, isdir
from os import getcwd
from pathlib import Path
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, TAXDB_URL, CDD_DATABASE_URL, TAXDB_SHA256
from celery_blast.dataset_refresh import DatasetRefreshSpec, refresh_dataset
from celery_blast.processes import ExternalCommandError, ExternalCommandTimeout, run_external_command


def build_startup_taxdb_download_command(taxdb_ftp_path, output_path):
    return ["wget", taxdb_ftp_path, "-q", "-O", output_path]


def build_startup_taxdb_extract_command(taxdb_archive_path, database_dir):
    return ["tar", "-zxvf", taxdb_archive_path, "-C", database_dir]


def run_startup_command(command, timeout):
    return run_external_command(command, timeout=timeout, shell=False, check=True).returncode


def startup_taxdb_refresh_spec():
    return DatasetRefreshSpec(
        name="taxdb",
        source_url=TAXDB_URL,
        public_root=Path(BLAST_DATABASE_DIR),
        required_files=("taxdb.btd", "taxdb.bti"),
        expected_sha256=TAXDB_SHA256,
        archive_name="taxdb.tar.gz",
        archive_type="tar.gz",
        expose_as="files",
        minimum_total_bytes=1,
        timeout=600,
    )


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

        check_domain_database_status()
    else:
        domain_database = DomainDatabase.objects.all()[0]
        if domain_database.domain_database_loaded == False:
            print("INFO: Conserved Domain Database is not loaded. CATHI setup has to be done by the user.")
            check_domain_database_status()
        elif domain_database.domain_database_loaded == True:
            if isdir(BLAST_DATABASE_DIR + 'CDD') == False:
                print("INFO: CDD database directory does not exist.")
                print("INFO: checking status of domain database ...")
                check_domain_database_status()

    if isfile(BLAST_DATABASE_DIR + 'taxdb.btd') and isfile(BLAST_DATABASE_DIR + 'taxdb.bti'):
        print("INFO:TAXONOMY DATABASE IS LOADED")

    else:
        print("INFO:NO TAXONOMY DATABASE")
        print("INFO:STARTING TO DOWNLOAD TAXONOMY DATABASE")

        try:
            current_working_directory = getcwd()  # /blast/reciprocal_blast
            print("INFO:TAXDB_URL: {}".format(TAXDB_URL))
            print("INFO:WORKING_DIRECTORY: {}".format(current_working_directory))

            metadata = refresh_dataset(startup_taxdb_refresh_spec())
            print("INFO:DONE DOWNLOADING TAXONOMY DATABASE version {}".format(metadata["version"]))

        except ExternalCommandTimeout as e:
            print("ERROR:TIMEOUT EXPIRED DURING DOWNLOAD OF TAXONOMY DATABASE: {}".format(e))
            print(
                "INFO:IF YOU HAVE NO STABLE INTERNET CONNECTION TRY TO RESTART THE WEBSERVER ONCE YOU HAVE A BETTER CONNECTION")
            print("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")
            print("ERROR:INSTALL THE TAXONOMY DATABASE MANUALLY into: {}".format("data/databases/"))
        except ExternalCommandError as e:
            print("ERROR:TAXONOMY DATABASE COMMAND RESULTED IN AN ERROR: {}".format(e))
            print("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")
            print("ERROR:INSTALL THE TAXONOMY DATABASE MANUALLY into: {}".format("data/databases/"))
