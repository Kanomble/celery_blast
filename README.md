# celery_blast
Symmetrical BLAST and target sequence search web interface with Django, Gunicorn, Ngninx, PostgreSQL, Celery, RabbitMQ, 
E-Direct, BLAST, Snakemake and Miniconda.

## Content
- [Installation](#installation)
- [Container network and configuration](#configuration_notes)
- [Project Setup](#project_setup)
- [BLAST Database creation](#blast_database)
  - [Notes on the download procedure](#download_process)
  - [Notes on the BLAST database formatting procedure](#makeblastdb)
- [Upload your own genome files](#genome_upload)
- [Result Dashboard](#result_dashboard)
- [Technical Details](#technical_details)
<a name="installation"></a>
## Installation
The application can get installed by submitting the `docker-compose up` or `docker compose up` command in a terminal window,
which points to the applications working directory (directory with `docker-compose.yml`). 
The docker client will pull remotely available images, the base image for this application, 
an image for the PostgreSQL database and finally an image for the RabbitMQ message broker.
The SymBLAST base image is pulled from this [DockerHub](https://hub.docker.com/repository/docker/kanomble/rec_blast_base_image).
All required software tools are loaded and installed automatically, there is no need for other third party software.

The base image has been build by using the Dockerfile of this repository.
Due to changes to some particular packages, the local build process may vary.
It is recommended to install the application with the remotely available images.

Next, docker creates seven containers named: `celery_blast_X_1` where `X` is a synonym for 
`nginx, worker, flower, web, postgres and rabbitmq`.

If you want to rebuild your docker images due to some (maybe fixed) error consider the cmd `docker-compose up --build` which will trigger a rebuild process (based on the context).
The web container will automatically try to restart if the startup fails, unless it is stopped manually (e.g. with Docker Desktop).

<a name="configuration_notes"></a>
### Notes on SymBLAST containers and possible configurations
SymBLAST is a server site tool, by starting the container network, your local computer will be used as a web
server. `Django` is the underlying web-framework and `gunicorn` serves as the WSGI HTTP Server. Both applications reside
in the SymBLAST base image. `Nginx` is used as a reverse proxy server, it directs client requests to `gunicorn`.
The long-running background tasks are managed by `rabbitmq` and `celery`, thus triggered processes are picked up by 
the message broker `rabbitmq` and passed to a queue, if a `celery-worker` is free, the process is executed. The task progress
is saved within the `postgresql` database within the `django_celery_results_taskresult` table, which enables task monitoring.
The `flower` container can be used to monitor the `celery-worker`. The reciprocal BLAST pipeline and the normal 
one-way BLAST pipelines are integrated into a Snakefile, which is used by the workflow management system `snakemake`.
Customization of Snakefiles enables user defined post-processing. In addtion, a `jupyter-notebook` container is 
integrated into the SymBLAST container network. Configuration is done within the `.env.prod` file. 
All important environment variables are defined within this file 
(e.g. the `DJANGO_ALLOWED_HOSTS` and the `SECRET_KEY` variables).

<a name="project_setup"></a>
## Project setup
To execute the integrated reciprocal BLAST pipeline of SymBLAST, certain data must be set up by the user. 
This includes the query sequences from a particular organism/genome file,
a forward BLAST database that will serve as the search space, a backward BLAST database, the scientific name of the 
organism from which the query sequences were obtained, and a project title. Additionally, the user can modify some BLAST
settings, such as the number of output sequences per query sequence (num_alignments) or the e-value cut-off. 
The BLAST databases can be selected from a special drop-down menu.

The forward BLAST database acts as a search space in which putative orthologous sequences can be located, 
while the backward BLAST database must contain the genome file from which the query sequences were obtained. 
Prior to saving a project into the database or executing the pipeline, the user-provided data undergoes validation to 
ensure that it meets the necessary criteria. If any of the validations fail,
accurate error messages are displayed within the relevant form fields to ensure a smooth pipeline execution.

The pipeline comprises the following steps:

1. Forward BLAST (default BLAST settings: e-value=0.001, word-size=3, threads=1, num_alignments=10000, max_hsps=500)
2. Backward BLAST preparation (extracting homologous target sequences of the forward BLAST)
3. Backward BLAST (BLAST search of the homologous target sequences against the genome of the query sequences)
4. Extracting Reciprocal Best Hits (RBHs, this is done via pandas merging tools)
5. Post-processing of RBHs (inference of taxonomic information, statistics, HTML and CSV tables, basic result plots)
6. Extraction of RBH-sequences separated by query sequences
7. Multiple sequence alignment of each set of RBHs with MAFFT
8. Phylogenetic inference of each set of RBHs with FastTree
9. Post-processing of the phylogenetic tree with ete3
10. CDD domain search of target sequences

Further SymBLAST post-processing procedures outside the scope of this pipeline involves:

1. Combining taxonomic information of the underlying database with the RBH result table
2. Building an interactive bokeh plot, that enables intuitive result interpretation
   1. Filter RBHs based on taxonomy, e-value, bitscore, sequence length and percent identity
   2. Download a selection of protein identifier
3. Refined phylogenetic inference with CDD domains of the RBHs
   1. Conducting RPS-BLAST with a specified set of RBHs
   2. Conducting a principal component analysis (PCA) based on the percent identity of the inferred domains with respect to the query sequence domains
   3. Building an interactive bokeh plot with the first two principal components and the taxonomic information within the RBH result table
   4. Refine the 

### Interactive Bokeh Plot
![Interactive Bokeh Plot](./celery_blast/static/images/example_bokeh_plot.png)

### Best practices for project settings
Use appropriate BLAST databases. If you want to search in more complete genomes, create a database that contains genome sequences
with a completeness level of `Chromosome` or `Complete Genome`. The `e-value` is more accurate for bigger databases, adjust the 
`e-value` according to your needs, this may have a huge effect on your inferred RBHs. Adjust the `num_alignments` parameter if
you work with huge databases, especially if you are working with more common sequences. 

The backward BLAST database should contain only one genome that corresponds to the taxonomic unit translated provided scientific name.


<a name="blast_database"></a>
## BLAST Databases
## BLAST database preparation
First, you need to tell SymBLAST to download the refseq assembly summary file from the
refseq [FTP](ftp://ftp.ncbi.nih.gov/genomes/refseq/) directory.
The application loads the summary file into a pandas dataframe, 
that is processed according to the user specifications. BLAST database user specifications are
the level of assembly completeness (e.g. 'Complete Genome', 'Chromosome', 'Contig' and 'Scaffold') 
and taxonomic information.

E.g. if you want to create a high quality database containing only species from the phylum `Cyanobacteriota`, you have to specify 
this during database creation. This can be done by typing `Cyanobacteriota` into the field "Scientific Names (sep. by ",")" and by
checking the assembly levels `Complete Genome` and `Chromosome`. 

If the user submits the form, a `BlastDatabase` model instance and a 
database directory, with a csv file containing the database table, is created.
The model is saved into the database, the database is not downloaded and formatted directly. 
The download and formatting procedure has to be started separately, which enables the user to validate database entries.
The download and format process progression is visualized on the database dashboard.
Available databases are shared between users.

Protein sequence files are downloaded from the NCBI FTP site and are passed to the `makeblastdb` command. Every 

<a name="download_process"></a>
## Download and formatting procedure of genome assemblies
When the user clicks the download button, a celery asynchronous task is initiated. 
This task consists of multiple subtasks, each responsible for performing specific steps in the database creation process. 
Celery allows for task monitoring through a logger and the ``TaskResult`` model of the [result backend](https://docs.celeryproject.org/en/stable/userguide/tasks.html#task-result-backends).

The task is triggered by a function called ``download_blast_databases_based_on_summary_file``
in the ``refseq_transactions/tasks.py`` file, which has three main steps. Firstly, it attempts to download and decompress 
all files stored in a table created during ``BlastDatabase`` model creation (as described in the above section). 
The celery task then formats the downloaded protein fasta files into BLAST databases and generates
an alias file similar to the one created with the ``blastdb_aliastool``. 
All three main steps are performed using the ``subprocess.Popen`` interface.

In the first step, the Unix commands ``wget`` and ``gzip`` are used to download and simultaneously 
decompress genome assemblies from the NCBI FTP server. This is accomplished with the following command: 
``wget -qO- ftp_path | gzip -d > assembly_file.faa``. 
If this command produces an error (returncode != 0), it is retried up to a maximum of ten times. 
If the command still fails after ten attempts, the assembly file is skipped, and the FTP path pointing to the file is written to a logfile.
Finally, in the next step, ``makeblastdb`` is used to format the downloaded and decompressed assembly files into BLAST databases.


<a name="makeblastdb"></a>
#### Notes on the `makeblastdb` program and BLAST database formatting
The ``makeblastdb`` module is a useful tool for building custom BLAST databases.
By default, it creates a database based on the input sequences. 
For example, if you submit the following command: ``makeblastdb -in .\prot_1_db.faa -dbtype prot -taxid 1140 -blastdb_version 5``,
a database named bw_prot_db.faa will be created. 
The ``-taxid`` parameter is used to assign the taxonomic node 1140 (corresponding to Synechococcus elongatus 7492) 
to all sequences in the bw_prot_db.faa fasta file.
If you need to format multiple fasta files using the ``makeblastdb`` program, there are two options available. 
First, you can provide multiple fasta files to the ``-in`` parameter. 
Alternatively, after formatting each fasta file individually, you can create a ``.pal`` database alias file that lists 
all existing databases. The ``blastdb_aliastool`` program can also be used to create this alias file for you.
If the sequences in the fasta files have the same taxonomic node, you can use the ``-taxid`` parameter to assign the taxid to the program. 
However, if they have different taxonomic nodes, you should use the ``-taxid_map`` parameter instead.

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

<a name="genome_upload"></a>
## Uploading genomes for BLAST database creation
The second option to obtain BLAST databases is to upload your own genomes. 
Currently, only genomes with protein sequences are supported. There are two different forms that can be used 
for uploading your own genome files.
1. The first form allows uploading a concatenated genome FASTA files with 
metadata file fields such as a taxmap file, which holds taxonomic information, an organism file, an assembly level file and
an assembly accession file. Most of these files are not mandatory. 
2. The other form allows uploading of multiple, single genome files together with their valid scientific organism names.

<a name="result_dashboard"></a>
## Result Dashboard
After successfully pipeline execution, results are displayed within a project details page.

<a name="snakemake"></a>
## Snakemake pipeline with celery
In order to enable reproducibility and an easy-to-use workflow execution, the workflow engine snakemake is used.
Snakemake associated snakefiles reside in a static directory `celery_blast/celery_blast/static/snakefiles`. 
They can be used outside this application, e.g. if the researcher needs to use additional settings or want to implement own post-processing procedures.
Different snakefiles are designed to execute the desired workflows. 
Currently, there are Snakefiles for the One-Way BLAST remote and local searches and the reciprocal BLAST analysis.
Execution of snakemake is wrapped in functions within the `tasks.py` files,
which are decorated with the celery `@shared_task` decorator. Those functions are queued up by rabbitmq and are processed by the celery worker.
Snakemake is executed with the `subprocess.Popen` interface which spawns a child process for every Snakemake workflow. 
Snakemake's working directory is the current project directory, e.g. if you executed the snakemake pipeline for a reciprocal BLAST
project, the `CWD` of snakemake is the path that points to this project. Default project paths are defined in the `celery_blast/settings.py` file.

Detailed information about Snakemake can be found on the [project documentation page](https://snakemake.readthedocs.io/en/stable/).

<a name="technical_details"></a>
## Technical Details
### Docker
All necessary software packages are deployed within docker container. Those containers are wrapped into a network
using docker-compose.

### Django, Gunicorn and Nginx
Django, Gunicorn, and Nginx are commonly used technologies for building and deploying web applications. 
In a Docker network, these technologies can be used together to create a scalable and efficient web application stack.
Django is a popular Python-based web framework used for building web applications. 
Gunicorn is a Python WSGI HTTP server that can serve Django applications.
Nginx is a high-performance web server that can act as a reverse proxy, load balancer, and serve static files.
When using Django, Gunicorn, and Nginx in a Docker network, each technology can be run in a separate Docker container. 
The Django application can be run using Gunicorn as a WSGI HTTP server, while Nginx can be used as a reverse proxy to direct traffic to the appropriate container.
This setup can be scaled horizontally by adding more containers to handle increased traffic. 
In addition, Docker's networking features can be used to ensure that traffic is properly directed to the appropriate container.
Overall, using Django, Gunicorn, and Nginx in a Docker network can create a highly scalable and efficient web application stack.

### POSTGRESQL database models
Django models reside inside project specific `models.py` files. Models are translated to database tables.
Documentation about the django.db.models package can be found [here](https://docs.djangoproject.com/en/2.2/topics/db/models/).
Relationships between models are managed with django model functions.
The `models.ForeignKey()` function is used for OneToMany / ManyToOne relations.
Additionally, there are the `models.OneToOneField()` and `models.ManyToManyField()` functions.
Relationships can get further described with `related_name` and `related_query_name` parameters, described in
[this](https://docs.djangoproject.com/en/2.2/ref/models/fields/#django.db.models.ForeignKey.related_query_name) Django documentation section.

Model managers reside inside project specific `managers.py` files.
Manager classes are responsible for the initial creation of the database models, such as `create_blast_project(fields=values ...)`.
Those functions can be used to trigger side effects during initialization of the database entry.
E.g. creation of blast project directories or specific setting files, such as the snakemake configuration file.

## TODO
- [ ] check if uploaded genomes consist of protein sequences
- [ ] refactor entrez query for one way blast remote searches
- [ ] deletion of taxonomic nodes
- [ ] interface for the CDD domain search
- [ ] separate download buttons for phylogenies and multiple sequence alignments within the external tools website
- [ ] database statistics progress bar
- [ ] duplicate entries in CDD result dataframe for bokeh plot - due to different assemblies - maybe too much info
- [ ] gunicorn --timeout setting
- [ ] section for each database statistic table
- [ ] update BLAST databases
- [ ] if no genome level is defined take all
- [ ] snakemake --unlock ? 
- [ ] refactor website structures
- [X] add a production environment
  - [X] follow this [guide](https://testdriven.io/blog/dockerizing-django-with-postgres-gunicorn-and-nginx/)
- [X] add a development environment
- [X] refactor one-way-remote BLAST bokeh plots
  - [X] download button
  - [X] annotations in plot
  - [X] add progress bar
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
- [X] fasttree, mafft timeout setting: 4000s not appropriate: 40.000s
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


## useful documentation:
- Interaction with NCBI (Entrez) via python [Biopython package](https://biopython.org/wiki/Documentation)
- [Celery Project Documentation](https://docs.celeryproject.org/en/stable/django/first-steps-with-django.html)
- Documentation for [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
- Documentation for [celery-progress](https://github.com/czue/celery-progress) - youtube tutorial [celery-progress](https://www.youtube.com/watch?v=BbPswIqn2VI)
- [BLAST DB](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) FTP server description
- [E-Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
