'''
check out the official django documentation for informations regarding this
python file - https://www.djangoproject.com/
'''

from django.contrib import admin

# Register your models here.
from .models import BlastProject, BlastDatabase, BlastSettings, AssemblyLevels

admin.site.register(BlastProject)
admin.site.register(BlastDatabase)
admin.site.register(BlastSettings)
admin.site.register(AssemblyLevels)