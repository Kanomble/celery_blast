# celery_blast
Reciprocal BLAST web-interface with Django, Celery, Flower, RabbitMQ, E-Direct and BLAST
## TODO
- [ ] correct timezone in the docker image
- [x] integrate functionality for Create Taxonomic Node File option in celery_blast project
    - [ ] think about multiple species_name inputs ...
- [ ] integrate blast_project_dashboard functionality
    - [x] TODO Database Models for dashboard functionality (base functionality)
    - [ ] add links to Detail View, Delete View and Execution View
- [ ] check out the BLAST Database software
- [ ] check out the .pal files from BLAST databases

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

Wrap database transactions inside `with transactions.atomic()` blocks. 

## blast_project dashboard


## useful documentation:
- Interaction with NCBI (Entrez) via python [Biopython package](https://biopython.org/wiki/Documentation)
- [Celery Project Documentation](https://docs.celeryproject.org/en/stable/django/first-steps-with-django.html)
- [BLAST DB](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) FTP server description
