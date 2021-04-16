from .models import BlastProject, BlastDatabase

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