# celery_blast
Reciprocal BLAST web-interface with Django, Celery, Flower, RabbitMQ, E-Direct and BLAST
## Installation
Installation can be done with `docker-compose up`. To activate the E-Direct tool do following:
```` Bash
docker exec -it celery_blast /bin/bash
#in the docker shell:
cd ../edirect && sh ./setup.sh
#answer with y
````
## TODO
- [ ] write documentation for added functions
- [ ] refactor the create_blastdatabase_table_and_directory function (too long)
- [X] refactor the refseq_transactions_dashboard in order to allow creation of a database directory with a csv table file, and tables for deletion download and details functions
    - [ ] add table for not downloaded databases with delete and download button (download button triggers snakemake)
    - [X] add table for downloaded databases with deletion button
    - [X] add table for downloaded databases with errors and with a deletion button
- [ ] correct timezone in the docker image
- [x] integrate functionality for Create Taxonomic Node File option in celery_blast project
    - [ ] think about multiple species_name inputs ...
- [ ] integrate blast_project_dashboard functionality
    - [x] TODO Database Models for dashboard functionality (base functionality)
    - [ ] add links to Detail View, Delete View and Execution View
- [ ] integrate blastdb_creation dashboard functionality
    - [ ] assembly download, `makeblastdb` execution and creation of `.pal` blast database alias files with snakemake
    - [ ] check out how to monitor snakemake execution
- [ ] integrate project_creation dashboard functionality
- [X] use [ajax](https://api.jquery.com/jquery.ajax/) asynchronously
    - [ ] use [ajax](https://api.jquery.com/jquery.ajax/) asynchronously to check celery tasks execution process
- [ ] use [celery-progress](https://github.com/czue/celery-progress) for monitoring the celery tasks execution process in the backend
- [X] check out the .pal files from BLAST databases

## TODO Database Models
- [ ] create models:
    - [X] BlastProject
    - [X] [BlastProjectManager](https://docs.djangoproject.com/en/2.2/ref/models/instances/)
    - [X] BlastDatabase
    - [ ] BlastDatabaseManager
    - [X] BlastSettings
    - [X] AssemblyLevels
- [ ] add validation
- [ ] write tests
- [ ] refactor models
- [ ] wrap database transactions inside `with transactions.atomic()` blocks

## database models
Django models reside inside project specific `models.py` files. Models are translated to database tables. 
Documentation about the django.db.models package can be found [here](https://docs.djangoproject.com/en/2.2/topics/db/models/).
Relationships between models are managed with django model functions.
The `models.ForeignKey()` function is used for OneToMany / ManyToOne relations. 
Additionally there are the `models.OneToOneField()` and `models.ManyToManyField()` functions. 
Relationships can get further described with `related_name` and `related_query_name` parameters, described in this
[this](https://docs.djangoproject.com/en/2.2/ref/models/fields/#django.db.models.ForeignKey.related_query_name) Django documentation section.

Model managers reside inside project specific `managers.py` files. 
Manager classes are responsible for specific creation functions, such as ``create_blast_project(fields=values ...)``.
Those functions can be used to trigger side effect during initialization of the database entry.
E.g. creation of blast project directories or filepath settings...

## BLAST Databases
### BLAST database preparation
Currently BLAST databases are downloaded from the BLAST FTP site provided by the NCBI. First the software downloads the refseq assembly summary file from the refseq [FTP](ftp://ftp.ncbi.nih.gov/genomes/refseq/) directory. This summary file inherits 226337 assembly entries. The application loads this summary file into a pandas dataframe, that gets processed. As a first step of database creation, the user has to define the level of assembly completeness (e.g. 'Complete Genome', 'Chromosome', 'Contig' and 'Scaffold') and optionally a taxonomic node file. Based on the `assembly_level` and (if the taxonomic node file is provided) `taxid` columns the summary file gets filtered. For example, the user could specify the assembly levels of the new database as `Complete Genome` and `Chromosome` and the `apes.taxids` file as taxonomy limitation. According to this setup, the summary file gets filtered by the provided taxids (which reside in the `apes.taxids` file) and the assembly levels, which results into an table with 6 entries (20.04.2021).

If the user submits the form, a `BlastDatabase` model instance and a database directory with a file, that contains the filtered table, is created. The model is saved into the database without an associated `TaskResult`, thus yet it is not downloaded and formatted. 

### Formatting procedure of downloaded genome assemblies
If the user press the download button, a celery asynchronous task is executed, which in turn executes a snakemake process via the `subprocess.Popen` interface.
The snakefile that gets executed is located in a static folder, but snakemake's working directory will be set to the `BlastDatabase` filepath.

With the `makeblastdb` module we can create custom BLAST databases. The `makeblastdb` cmd is executed inside a snakemake script, that is executed as a celery `@shared task`.
Per default it will create a database for the input sequences, e.g. if you submit following cmd: 
`makeblastdb -in .\prot_1_db.faa -dbtype prot -taxid 1140 -blastdb_version 5` you will create a database for the bw_prot_db.faa.
The `-taxid` parameter is used to assign the taxonomic node 1140 (*Synechococcus elongatus* 7492) to all sequences that reside in the `bw_prot_db.faa` fasta file.
If you have mutliple fasta files that needs to get formatted with the `makeblastdb` program, there are two options. First, you can path multiple fasta files to the `-in` parameter. Secondly, after formatting each fasta individually, you can create a `.pal` database alias file that lists all existing databases or you can use the `blastdb_aliastool` program to create this alias file for you. If they have the same taxonomic node you can parse the taxid with the `-taxid` parameter to the program, if not use the `-taxid_map` parameter. 

```` Bash
makeblastdb -in .\prot_1_db.faa -dbtype prot -taxid 1140 -blastdb_version 5
makeblastdb -in .\prot_2_db.faa -dbtype prot -taxid 1844971 -blastdb_version 5
blastdb_aliastool -dblist 'prot_1_db.faa prot_2_db.faa' -dbtype prot -title combined_db -out combined_db
blastp -query .\test.faa -db combined_db -out blast_out.table
````
Example of the `combined_db.pal` file:
````Text
#
# Alias file created 04/17/2021 12:50:29
#
TITLE combined_db
DBLIST "prot_1_db.faa" "prot_2_db.faa"
````
The `.pal` file combines different formatted BLAST databases so that they can be used like one combined database. 
This is useful for databases with duplicate sequences, they normally have an identifier (accession number) that starts with `WP`.

## SNAKEMAKE tasks with celery
In order to allow reproducability and allow an easy workflow understanding, the workflow engine snakemake is used. 
Snakemake associated snakefiles reside in a static directory `celery_blast/celery_blast/static/`. 
Different snakefiles are designed to execute the desired workflow. Execution of snakemake is wrapped in functions of the `tasks.py` files,
which are decorated with the celery `@shared_task` decorator. Those function use the `subprocess.Popen` interface to spawn the snakemake process.
During execution the underlying database (e.g. BlastDatabase or BlastProject) model OneToOne field gets updated with the appropriate `TaskResult` model.
This allows interaction with the associated celery task and can be used for displaying the progress of the task. 
Furthermore, snakemake is executed with the `--wms-monitor` parameter, that enables snakemake communication with [Panoptes](https://github.com/panoptes-organization/monitor-schema). In addition [Flower](https://flower.readthedocs.io/en/latest/) can be used to monitor the celery tasks.
### TODO snakemake
- [X] design a snakefile for downloading blast databases
- [ ] design a snakefile for reciprocal BLAST analysis
- [ ] messages during tasks execution to [celery-progress](https://github.com/czue/celery-progress)

## POSTGRESQL database transactions

## blast_project dashboard


## useful documentation:
- Interaction with NCBI (Entrez) via python [Biopython package](https://biopython.org/wiki/Documentation)
- [Celery Project Documentation](https://docs.celeryproject.org/en/stable/django/first-steps-with-django.html)
- Documentation for [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
- Documentation for [celery-progress](https://github.com/czue/celery-progress)
- [BLAST DB](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) FTP server description
