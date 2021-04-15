# celery_blast
Reciprocal BLAST web-interface with Django, Celery, Flower, RabbitMQ, E-Direct and BLAST
## TODO
- [x] integrate functionality for Create Taxonomic Node File option in celery_blast project
    - [ ] think about multiple species_name inputs ...
- [ ] create models:
    - BlastProject
    - BlastDatabase
    - BlastSettings

## database models

````python
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult

class BlastProject(models.Model):
    project_title = models.CharField(max_length=200, blank=False, unique=True)
    search_strategy = models.CharField(max_length=20, choices=[('blastp', 'blastp'), ('blastn', 'blastn')], default='blastp')
    project_username = models.ForeignKey(User,on_delete=models.CASCADE)
    project_database = models.ForeignKey(BlastDatabase,on_delete=models.CASCADE)

class BlastDatabase(models.Model):
    database_name = models.CharField(max_length=200, blank=False, unique=True)
    database_download_and_format_task = models.OneToOneField(TaskResult,on_delete=models.CASCADE,blank=True,null=True)
    database_description = models.CharField(max_length=200,unique=True)
    attached_taxonomic_node_file = models.CharField(max_length=300,blank=True,null=True)
    path_to_file = models.CharField(max_length=300,blank=True,null=True)
    #use the assembly_levels.SQL script for uploading the four existing assembly levels into the database
    assembly_levels = models.ManyToManyField(to=AssemblyLevels)
    assembly_entries = models.IntegerField()
    timestamp = models.DateTimeField(auto_now=True)

...
````