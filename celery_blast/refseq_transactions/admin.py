from django.contrib import admin

# Register your models here.
from .models import BlastDatabase, AssemblyLevels

admin.site.register(BlastDatabase)
admin.site.register(AssemblyLevels)
