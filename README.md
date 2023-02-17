# celery_blast
Reciprocal BLAST web-interface with Django, Celery, Flower, RabbitMQ, E-Direct, BLAST, Snakemake and Panoptes.
## Content

## Installation
The application can get installed by submitting the `docker-compose up` command in a terminal window which points to the applications working directory (directory with `docker-compose.yml`). The docker client will pull remotely available images, the base image for this application, an image for the PostgreSQL database and finally an image for the RabbitMQ message broker.
All images are pulled from this [DockerHub](https://hub.docker.com/repository/docker/kanomble/rec_blast_base_image).

The base image has been build by using the Dockerfile of this repository, a local build process
for the base image can get triggered by using the `docker-compose up -f installation_with_dockerfile_docker-compose.yml` command. 
All required software tools are loaded and installed automatically into this third image.
However, due to changes to some particular packages, the local build process can run into errors, especially the installation of `snakemake` will take some time and consume a lot of memory. 
It is recommended to install the application with the remotely available images.

Next, docker creates six containers named: `celery_blast_X_1` where `X` is a synonym for `panoptes, worker, flower, web, postgres and rabbitmq`.

If you run into any error during `docker-compose up` or if you recreate the container you need to activate the E-Direct tool. This can be achieved by submitting following command inside the docker container:
```` Bash
docker exec -it celery_blast /bin/bash
#in the docker shell:
cd ../edirect && sh ./setup.sh
#answer with y
````

If you want to rebuild your docker images due to some (maybe fixed) error consider the cmd `docker-compose up --build` which will trigger a rebuild process (based on the context).
The web container will automatically try to restart if the startup fails, unless it is stopped manually (e.g. with Docker Desktop).

## TODO
- [ ] refactor website structures
- [ ] add a production environment
  - [ ] follow this [guide](https://testdriven.io/blog/dockerizing-django-with-postgres-gunicorn-and-nginx/)
- [ ] add a development environment
- [ ] refactor one-way-remote BLAST bokeh plots
  - [ ] download button
  - [ ] annotations in plot
  - [ ] add progress bar
- [ ] add query sequence information to one way blast detail page - maybe also to table
- [ ] fix plotting: /blast/reciprocal_blast/media/blast_projects/2/.snakemake/scripts/tmp0ybykovj.build_folders_with_hit_info_for_each_qseqid.py:89: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).
  fig, ax = plt.subplots(2, 2)
- [ ] refactor correct deletion of static files - if database gets deleted, also delete project dirs associated to this database
- [ ] exception is thrown if there is no result in the esearch output
- [ ] refactor blast_tables_to_orthologous_table.py with new entrez function of the database statistics ncbi_transaction script
- [ ] refactor query_sequences_to_html_table.py - what happens if there are no information on NCBI available?
  - [ ] test with uploaded genomes and sequences
- [ ] failure BLAST database task
- [X] delete .gitkeep in postgres folder and fix wait-for script line ending
  - [X] wait-for script checking - LF / CLRF
  - [ ] alternatively write own .sh script file that checks if database/rabbitmq is loaded
- [ ] allow stopping database downloading and formatting tasks, allow deletion of paused/stopped database download/format procedures
- [ ] include taxonomic information in BLAST database tables
- [ ] esearch output into subfolders for each user
- [ ] there are some bugs in the EntrezSearch - e.g. if you search for rhino in pubmed
- [ ] what happens if a user deletes a database that is still connected to some projects?
- [ ] remote BLAST without entrez query: results in an error: Error: [blastp] internal_error: Message ID#54 Error: Failed to process the Entrez query: Only organism entrez queries are supported
- [ ] snakemake incomplete flag - restart function for snakemake
- [ ] limit for query sequences in one-way-blast
- [ ] add genbank database option
- [ ] upload genome: if \n is in any uploaded txt file it will count as a value for insertion
- [ ] display warning if backward organism not in the forward database as there is no controlling step - cause hits against the identical protein are not considered
- [ ] uploaded databases might have problems with taxonomic nodes - especially if the user selects different databases for the forward and backward blast
  - [ ] if there is no taxonomic node available (which can be the case for some organisms) it is not possible to upload a taxmap file ...
  - [ ] provide a tool for writing taxmap files of combined genomes 
  - [ ] forward and backward BLAST databases should both contain the backward BLAST taxonomic node, otherwise you wont see your 100% results in the backward BLAST!
- [ ] with uploaded DNA / RNA sequences it is currently not possible to perform the reciprocal BLAST analysis
    - [ ] multiple sequence alignment and phylogenetic tree reconstruction is also not supported with DNA/RNA seqs
    - [ ] refactor the reciprocal BLAST creation view - model supports blastn command (search strategy)
    - [ ] integrate the possible execution of blastn also in snakemake
- [ ] update taxdb option!
- [ ] input fasta query file headers have to be separated by a space, not by a pipe symbol (|) or others
  - [ ] reformat input sequences during upload
- [ ] exclude not downloaded and formatted assemblies from summary table
- [ ] write documentation for added functions
- [ ] blastn one way searches can't display query sequence information (of DNA sequences) received by biopython, biopython uses the protein db per default which causes errors if gene ids are provided
- [X] fix javascript download button for the protein database in the entrez search
- [X] add taxonomic information processing to database creation or database statistics calculation
- [X] refactor the external project information dashboard
  - [X] replace MSA and Phylogeny buttons with the number of RBH's
- [X] pipeline log files
  - [X] view logfile content on pipeline dashboard
- [X] refactor protein identifier field in one way BLASTs
- [X] fix error handling during one way BLAST project setup
- [X] add automate redirection to one way BLAST page
- [X] refactor one way blast rules
  - [X] blast_results_to_plots_and_html_table
  - [X] change altair to bokeh plots
- [X] refactor html result tables
  - [x] reciprocal results -> index column
  - [X] database statistics -> decimal place
- [X] refactor one-way remote BLAST pipeline
- [X] refactor logging in exceptions
- [X] table column names in entrez esearch dashboard
- [X] integrate MAFFT and FastTree in the corresponding snakefiles
- [X] define global timeout variable --> use celery timeouts --> soft and hard timeouts
- [X] correct ajax requests if it results into an error
- [X] implement custom snakemake logfile that lists all of the executed functions
- [X] if no reciprocal hits are available for at least one gene the snakemake workflow will result into an error - rule extract_sequences
- [X] reciprocal_result.csv should also contain assembly accession for each hit and genus, family .. informations
  - [X] genus, family and other taxonomic information are integrated
  - [X] assembly accession id
- [X] add correct docker timezone (currently it is Europe/Berlin but still 2h difference)
- [X] add pident ssequence and other output options to the blast results
- [X] check genome upload option | all options - begin with the view function
- [x] download taxdb during build process of the docker image
- [x] check if backward organism is in database
- [X] refactor the create_blastdatabase_table_and_directory function (too long)
- [x] extract subject sequences from database (with blastdbcmd and orthologous sequence id list)
    - [x] integrate docker container for mafft and fasttree 
    - [x] perform msa with orthologous subject sequences
    - [x] build ml or neighbour joining trees from all msa's
  - [X] check if query sequences are in backward database
- [X] installation still requires the `assembly_levels.sql` SQL-Script which inserts the four assembly levels, search for automatic insertions by installation
- [X] add more options to BlastSettings - Alter BlastSettings model and forms
- [x] integrate functionality for Create Taxonomic Node File option in celery_blast project
    - [X] think about multiple species_name inputs ... --> not possible
- [X] add configuration environment variables for SNAKEMAKE - settings.py
- [X] use a config file for all configuration options, e.g. the panoptes - settings.py
- [X] reciprocal BLAST results, sequence in output and string of taxonomy; EXAMPLE: https://docs.google.com/spreadsheets/d/1EBwEp-C0ocCUBx3zVaCtxrMlTKGTZHx19lPq8q7W_bE/edit#gid=649158632
- [X] error handling if ftp_path does not exist e.g.: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/893/775/GCF_000893775.1_ViralProj70005/GCF_000893775.1_ViralProj70005_protein.faa.gz
- [X] refactor the refseq_transactions_dashboard in order to allow creation of a database directory with a csv table file, and tables for deletion download and details functions
    - [X] add table for not downloaded databases with delete and download button (download button triggers snakemake)
    - [X] add table for downloaded databases with deletion button
    - [X] add table for downloaded databases with errors and with a deletion button
- [x] correct timezone in the docker image
- [x] integrate blast_project_dashboard functionality
    - [x] database models for dashboard functionality (base functionality)
    - [x] add links to Detail View, Delete View and Execution View
- [X] integrate blastdb_creation dashboard functionality
    - [X] assembly download, `makeblastdb` execution and creation of `.pal` blast database alias files with celery
    - [X] monitoring with ajax and celery progress
- [X] integrate project_creation dashboard functionality
- [X] use [ajax](https://api.jquery.com/jquery.ajax/) asynchronously
    - [X] use [ajax](https://api.jquery.com/jquery.ajax/) asynchronously to check celery tasks execution process
- [X] use [celery-progress](https://github.com/czue/celery-progress) for monitoring the celery tasks execution process in the backend
- [X] check out the .pal files from BLAST databases

## TODO Database Models
- [X] add validation
- [X] write tests
- [X] wrap database transactions inside `with transactions.atomic()` blocks
- [X] create models:
    - [X] BlastProject
        - [X] add backward BlastDatabase
    - [X] [BlastProjectManager](https://docs.djangoproject.com/en/2.2/ref/models/instances/)
    - [X] BlastDatabase
    - [X] BlastDatabaseManager
    - [X] BlastSettings
    - [X] AssemblyLevels
    - [X] UploadedGenomes
    - [X] ExternalTools --> Connected to BlastProjects
    - [X] QuerySequences --> Connected to ExternalTools

## Project setup
In order to execute the integrated reciprocal BLAST pipeline (the SymBLAST core pipeline), the user has to
set up some essential data: the query sequences from one particular organism/genome file, a forward BLAST database, a
backward BLAST database, the scientific name of the organism from which the query sequences were obtained and a project title. 
In addition, the user may change some BLAST settings, e.g. the number of output sequences per query sequence (`num_alignments`),
or the `e-value` cut-off. The BLAST databases can be selected via a special drop-down menu.
The Forward BLAST database serves as a search space where putative orthologous sequences can be found.
The backward BLAST database has to contain the genome file from which the query sequences were obtained.
## BLAST Databases
## BLAST database preparation
Protein sequence files are downloaded from the NCBI FTP site and are passed to the `makeblastdb` command.
First the software downloads the refseq assembly summary file from the
refseq [FTP](ftp://ftp.ncbi.nih.gov/genomes/refseq/) directory.
The application loads the summary file into a pandas dataframe, 
that is processed according to the user specifications. BLAST database user specifications are
the level of assembly completeness (e.g. 'Complete Genome', 'Chromosome', 'Contig' and 'Scaffold') 
and taxonomic information.
For example, the user could specify the assembly levels 
of the new database as `Complete Genome` and `Chromosome` and the `apes.taxids`
file as basis for taxonomic limitation. 
According to this setup, the summary file gets filtered by the provided 
taxids (which reside in the `apes.taxids` file) and the assembly levels,
which results into a table with 6 entries (20.04.2021). Files for the taxonomic limitations can get uploaded by the user
or directly created with the SymBLAST tool. 

If the user submits the form, a `BlastDatabase` model instance and a 
database directory with a file, that contains the filtered table, is created.
The model is saved into the database, the database is not downloaded and formatted directly. 
The download and formatting procedure has to be started separately, which enables the user to validate database entries.
The download and format process progression is visualized on the database dashboard.
Available databases are shared between users.

## Download and formatting procedure of genome assemblies
If the user presses the download button, a celery asynchronous task is executed. 
This task is composed of multiple subtasks, that perform the relevant database creation steps.
Celery allows task monitoring via a logger and with the `TaskResult` model of the [result backend](https://docs.celeryproject.org/en/stable/userguide/tasks.html#task-result-backends). 
The name of the function and celery task, which is triggered by the button is `download_blast_databases_based_on_summary_file` resides in the `refseq_transactions/tasks.py` file. 
This function has to process three main steps, first it tries to download and decompress all files that 
reside in a table that has been created during `BlastDatabase` model creation (s. above section). 
The celery task formats the downloaded protein fasta files to BLAST databases and writes an alias file similar to the file that is created with the `blastdb_aliastool`. 
In all of these three main steps processes via the `subprocess.Popen` interface are executed. 
As a first step, the unix command `wget` and `gzip` are used to download and simultaneous decompress genome assemblies from the NCBI FTP server.
This is done with following command: `wget -qO- ftp_path | gzip -d > assembly_file.faa`. If this command results into any error (returncode != 0) it is retried with a maximum of ten tries. 
If there is no positive returncode after ten tries, the assembly file is skipped and the FTP-path that points to the file is written into a logfile.
In the next step `makeblastdb` is used to format the downloaded and decompressed assembly files.

#### Notes on the `makeblastdb` program and BLAST database formatting
With the `makeblastdb` module custom BLAST databases are builded. Per default it will create a database for the input sequences, e.g. if you submit following cmd:
`makeblastdb -in .\prot_1_db.faa -dbtype prot -taxid 1140 -blastdb_version 5` you will create a database for the `bw_prot_db.faa`.
The `-taxid` parameter is used to assign the taxonomic node 1140 (*Synechococcus elongatus* 7492) to all sequences that reside in the `bw_prot_db.faa` fasta file.
If you have mutliple fasta files that needs to get formatted with the `makeblastdb` program, there are two options. First, you can pass multiple fasta files to the `-in` parameter. 
Secondly, after formatting each fasta individually, you can create a `.pal` database alias file that lists all existing databases, 
or you can use the `blastdb_aliastool` program to create this alias file for you. If they have the same taxonomic node you can pass the taxid with the `-taxid` parameter to the program, if not use the `-taxid_map` parameter.

```` Bash
makeblastdb -in .\prot_1_db.faa -dbtype prot -taxid 1140 -blastdb_version 5
makeblastdb -in .\prot_2_db.faa -dbtype prot -taxid 1844971 -blastdb_version 5
blastdb_aliastool -dblist 'prot_1_db.faa prot_2_db.faa' -dbtype prot -title combined_db -out combined_db
blastp -query .\test.faa -db combined_db -out blast_out.table -outfmt "6 qseqid sseqid evalue bitscore qgi sgi sacc pident nident mismatch gaps qcovhsp staxids sscinames scomnames sskingdoms  stitle"
````
BLASTP `-outfmt 6` output formats are described [here](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6).

The third step creates an alias file for all BLAST databases that have been formatted in the previous step.
The alias file is the file, that is created by the `blastdb_aliastool`.

Example of the `combined_db.pal` file:
````Text
#
# Alias file created 04/17/2021 12:50:29
#
TITLE combined_db
DBLIST "prot_1_db.faa" "prot_2_db.faa"
````
The `.pal` file combines different formatted BLAST databases so that they can be used as one combined database. 
This is useful for databases with duplicate sequences, they normally have an identifier (accession number) that starts with `WP`.

During execution the underlying database (e.g. `BlastDatabase` or `BlastProject`) model OneToOne field gets updated with the appropriate celery `TaskResult` model.
This allows interaction with the associated celery task and can be used for displaying the progress of the task.

## Uploading genomes for BLAST database creation
The second option to obtain BLAST databases is to upload your own genomes. 
Currently, only genomes with protein sequences are supported. There are two different forms that can be used 
for the uploading your own genome files. The first form allows uploading a concatenated genome FASTA files with 
metadata file fields such as a taxmap file, which holds taxonomic information, an organism file, an assembly level file and
an assembly accession file. Most of these files are not mandatory, except for the taxmap_file when no taxonomic node is 
is present. The other form allows uploading of single genome files together with their valid scientific organism names.

## SNAKEMAKE tasks with celery
In order to enable reproducability and develop and easy-to-use workflow execution, the workflow engine snakemake is used.
Snakemake associated snakefiles reside in a static directory `celery_blast/celery_blast/static/`. 
In addition they can be used outside this application, e.g. if the researcher needs to use additional settings or want to implement own post-processing procedures.
Different snakefiles are designed to execute the desired workflows. Currently, there are Snakefiles for the One-Way BLAST remote and local searches and the reciprocal BLAST analysis.
Execution of snakemake is wrapped in functions within the `tasks.py` files,
which are decorated with the celery `@shared_task` decorator. Those functions are queued up by rabbitmq and are processed by the celery worker.
Snakemake is executed with the `subprocess.Popen` interface which spawns a child process for every Snakemake workflow.

Furthermore, snakemake is executed with the `--wms-monitor` parameter, that enables snakemake communication with [Panoptes](https://github.com/panoptes-organization/monitor-schema). In addition [Flower](https://flower.readthedocs.io/en/latest/) can be used to monitor the celery tasks.

## POSTGRESQL database transactions
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
E.g. creation of blast project directories or file settings...

## blast_project dashboard

## useful documentation:
- Interaction with NCBI (Entrez) via python [Biopython package](https://biopython.org/wiki/Documentation)
- [Celery Project Documentation](https://docs.celeryproject.org/en/stable/django/first-steps-with-django.html)
- Documentation for [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
- Documentation for [celery-progress](https://github.com/czue/celery-progress) - youtube tutorial [celery-progress](https://www.youtube.com/watch?v=BbPswIqn2VI)
- [BLAST DB](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) FTP server description
- [E-Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
