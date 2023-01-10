'''
check out the official django documentation for informations regarding this
python file - https://www.djangoproject.com/
'''

from django.contrib import admin

# Register your models here.
from .models import BlastProject, BlastSettings

admin.site.register(BlastProject)
admin.site.register(BlastSettings)
