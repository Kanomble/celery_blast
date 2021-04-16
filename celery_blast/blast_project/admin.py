from django.contrib import admin

# Register your models here.
from .models import BlastProject, BlastDatabase, BlastSettings, AssemblyLevels

admin.site.register(BlastProject)
admin.site.register(BlastDatabase)
admin.site.register(BlastSettings)
admin.site.register(AssemblyLevels)