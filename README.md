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
Django models that are translated to database tables. Documentation about the django.db.models package can be found [here](https://docs.djangoproject.com/en/2.2/topics/db/models/).

````python
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult

class BlastProject(models.Model):
    BLAST_SEARCH_PROGRAMS = [('blastp', 'blastp'), ('blastn', 'blastn')]
    
    #attribute fields
    project_title = models.CharField(
                verbose_name="project title", 
                max_length=200, blank=False, unique=True)
    search_strategy = models.CharField(
                verbose_name="BLAST program that is used inside snakemake execution",
                max_length=20, 
                choices=BLAST_SEARCH_PROGRAMS, 
                default='blastp')
    timestamp = models.DateTimeField(auto_now=True)
    
    #relationships
    project_user = models.ForeignKey(
                User,
                on_delete=models.CASCADE,
                verbose_name="user who created the project")
    
    project_forward_settings = models.OneToOneField(
                BlastSettings,
                on_delete=models.CASCADE,
                verbose_name="settings for the forward BLAST execution")
    project_backward_settings = models.OneToOneField(
                BlastSettings,
                on_delete=models.CASCADE,
                verbose_name="settings for the backward BLAST execution")   
                       
    #forward database and backward database
    #each project can have ONE forward and ONE backward BlastDatabase
    # --> two foreign keys
    project_forward_database = models.ForeignKey(
                BlastDatabase,
                on_delete=models.CASCADE,
                verbose_name="associated forward BLAST database")
    
    project_backward_database = models.ForeignKey(
                BlastDatabase,
                on_delete=models.CASCADE,
                verbose_name="associated backward BLAST database")
    #one to one relationship
    project_execution_snakemake_task = models.OneToOneField(
                TaskResult,
                on_delete=models.SET_NULL,
                blank=True,null=True,
                verbose_name="django_celery_results taskresult model for this project")
    
    #overwritten functions
    def __str__(self):
        return "Reciprocal BLAST Project, created {} by {} with fw db {} and bw db {}".format(
                self.timestamp,self.project_user.name,
                self.project_forward_database.database_name,
                self.project_backward_database.database_name)

    #functions
    def get_project_username(self):
        return self.project_user.name
    
    def get_project_useremail(self):
        return self.project_user.email
    
    def get_project_dir(self):
        return 'media/blast_projects/' + str(self.id)
    
    def if_executed_return_associated_taskresult_model(self):
        #executed
        if self.project_execution_snakemake_task != None:
            return self.project_execution_snakemake_task
        #not executed
        else:
            return None


class BlastDatabase(models.Model):
    #attribute fields
    database_name = models.CharField(
                max_length=200, 
                blank=False, unique=True,
                verbose_name="database name")
    database_description = models.CharField(
                max_length=200,
                verbose_name="short description of database purpose")
    
    assembly_entries = models.IntegerField(
                verbose_name="number of assembly entries that should get downloaded")
    timestamp = models.DateTimeField(
                auto_now=True,
                verbose_name="date of database creation")    

    #nullable fields
    #possibility to add a taxonomic file
    attached_taxonomic_node_file = models.CharField(
                max_length=300,
                blank=True,null=True,
                verbose_name="associated taxonomic file, which was used to limit assembly entries in db creation by taxids")
    path_to_database_file = models.CharField(
                max_length=300,
                blank=True,null=True,
                verbose_name="after makeblastdb task has finished this field is set automatically with the path to the BLAST database")
    
    #relationships
    database_download_and_format_task = models.OneToOneField(TaskResult,on_delete=models.CASCADE,blank=True,null=True)
    #use the assembly_levels.SQL script for uploading the four existing assembly levels into the database
    assembly_levels = models.ManyToManyField(to=AssemblyLevels)

class BlastSettings(models.Model):
    pass
...
````

## useful documentation:
- Interaction with NCBI (Entrez) via python [Biopython package](https://biopython.org/wiki/Documentation)
- [Celery Project Documentation](https://docs.celeryproject.org/en/stable/django/first-steps-with-django.html)