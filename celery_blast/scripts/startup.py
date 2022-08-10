from blast_project.models import AssemblyLevels
from os.path import isfile
from os import remove, getcwd
import subprocess
import psutil

def run():
    if len(AssemblyLevels.objects.all()) != 4:
        print("INFO:Inserting Assembly Levels")
        contig = AssemblyLevels(assembly_level='Contig')
        chromosome = AssemblyLevels(assembly_level='Chromosome')
        complete = AssemblyLevels(assembly_level='Complete Genome')
        scaffold = AssemblyLevels(assembly_level='Scaffold')
        contig.save()
        chromosome.save()
        complete.save()
        scaffold.save()
    if isfile('media/databases/taxdb.btd') and isfile('media/databases/taxdb.bti'):
        print("INFO:TAXONOMY DATABASE IS LOADED")
    else:
        print("INFO:NO TAXONOMY DATABASE")
        if isfile("media/databases/taxdb.tar.gz"):
            remove("media/databases/taxdb.tar.gz")
        if isfile("media/databases/taxdb.tar"):
            remove("media/databases/taxdb.tar")
        print("INFO:STARTING TO DOWNLOAD TAXONOMY DATABASE")

        try:
            taxdb_ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"
            current_working_directory = getcwd()  # /blast/reciprocal_blast
            path_to_taxdb_location = current_working_directory + '/media/databases/'
            path_to_taxdb_location = path_to_taxdb_location + 'taxdb.tar.gz'

            proc = subprocess.Popen(["wget", taxdb_ftp_path, "-q", "-O", path_to_taxdb_location], shell=False)
            returncode = proc.wait(timeout=600)
            if returncode != 0:
                raise subprocess.SubprocessError
            print("INFO:EXTRACTING TAXONOMY DB")

            proc = subprocess.Popen(["tar", "-zxvf", path_to_taxdb_location, "-C", "/blast/reciprocal_blast/media/databases/"], shell=False)
            returncode = proc.wait(timeout=600)
            if returncode != 0:
                raise subprocess.SubprocessError
            print("INFO:DONE DOWNLOADING TAXONOMY DATABASE")
        except subprocess.TimeoutExpired as e:
            print("ERROR:TIMEOUT EXPIRED DURING DOWNLOAD OF TAXONOMY DATABASE: {}".format(e))
            print(
                "INFO:IF YOU HAVE NO STABLE INTERNET CONNECTION TRY TO RESTART THE WEBSERVER ONCE YOU HAVE A BETTER CONNECTION")
            print("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")

            pid = proc.pid
            parent = psutil.Process(pid)

            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()

        except subprocess.SubprocessError as e:
            print("ERROR:WGET RESULTED IN AN ERROR: {}".format(e))
            print("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")

            pid = proc.pid
            parent = psutil.Process(pid)

            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()

