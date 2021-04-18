# celery_blast
Reciprocal BLAST web-interface with Django, Celery, Flower, RabbitMQ, E-Direct and BLAST
## TODO
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
- [X] check out the .pal files from BLAST databases

## TODO Database Models
- [X] create models:
    - BlastProject
    - [BlastProjectManager](https://docs.djangoproject.com/en/2.2/ref/models/instances/)
    - BlastDatabase
    - BlastSettings
    - AssemblyLevels
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
With the `makeblastdb` module we can create custom BLAST databases. The `makeblastdb` cmd is executed inside a snakemake script, that is executed as a celery `@shared task`.
Per default it will create a database for the input sequences, e.g. if you submit following cmd: 
`makeblastdb -in .\prot_1_db.faa -dbtype prot -taxid 1140 -blastdb_version 5` you will create a database for the bw_prot_db.faa.
The `-taxid` parameter is used to assign the taxonomic node 1140 (*Synechococcus elongatus* 7492) to all sequences that reside in the `bw_prot_db.faa` fasta file.
If you have mutliple fasta files that needs to get formatted with the `makeblastdb` program, there are two options. First, you can path multiple fasta files to the `-in` parameter. Secondly, after formatting each fasta individually, you can create a `.pal` database alias file that lists all existing databases or you can use the `blastdb_aliastool` program to create this alias file for you. If they have the same taxonomic node you can parse the taxid with the `-taxid` parameter to the program, if not use the `-taxid_map` parameter. 

```` Shell
makeblastdb -in .\prot_1_db.faa -dbtype prot -taxid 1140 -blastdb_version 5
makeblastdb -in .\prot_2_db.faa -dbtype prot -taxid 1844971 -blastdb_version 5
blastdb_aliastool -dblist 'prot_1_db.faa prot_2_db.faa' -dbtype prot -title combined_db -out combined_db
blastp -query .\test.faa -db combined_db -out blast_out.table
````

## POSTGRESQL database transactions

## blast_project dashboard


## useful documentation:
- Interaction with NCBI (Entrez) via python [Biopython package](https://biopython.org/wiki/Documentation)
- [Celery Project Documentation](https://docs.celeryproject.org/en/stable/django/first-steps-with-django.html)
- Documentation for [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
- [BLAST DB](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) FTP server description
