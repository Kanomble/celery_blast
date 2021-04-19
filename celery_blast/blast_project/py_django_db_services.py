from .models import BlastProject, BlastDatabase, AssemblyLevels
from django.db import IntegrityError

#following functions are utilized in the dashboard_view
''' get_users_blast_projects
    
    returns a query-set of blast_projects from the logged in user
    
    :param userid (request.user.id)
        :type int
    :returns django.db.models.query.QuerySet of BlastProjects
'''
def get_users_blast_projects(userid):
    return BlastProject.objects.get_blast_projects_by_userid(userid)

''' get_all_blast_databases

:returns all blastdatabases as a query-set
'''
def get_all_blast_databases():
    return BlastDatabase.objects.all()



def create_and_save_refseq_database_model(database_name,database_description,assembly_levels,assembly_entries,attached_taxonomic_file=None):
    try:

        #create model refseq genome objects (s. models.py file)
        #path_to_database_file = 'media/' + 'databases/' + 'refseq_databases/' + database_description.replace(' ','_').upper() + '.database.faa'
        if attached_taxonomic_file != None:
            new_refseq_genome = BlastDatabase.objects.create(
                database_name=database_name,
                database_description=database_description,
                assembly_entries=assembly_entries,
                attached_taxonomic_node_file=attached_taxonomic_file)
        else:
            new_refseq_genome = BlastDatabase.objects.create(
                database_name=database_name,
                database_description=database_description,
                assembly_entries=assembly_entries)


        #get all associated assembly levels (max. 4)
        assembly_levels_models = AssemblyLevels.objects.filter(assembly_level__in=assembly_levels)

        for assembly_level in assembly_levels_models:
            new_refseq_genome.assembly_levels.add(assembly_level)

        new_refseq_genome.path_to_database_file = 'media/databases/' + str(new_refseq_genome.id)
        new_refseq_genome.save()
        return new_refseq_genome
    except Exception as e:
        raise IntegrityError('couldnt save refseq genome model into database with exception : {}'.format(e))