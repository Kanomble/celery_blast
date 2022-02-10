#from __future__ import absolute_import, unicode_literals
import os
from celery import Celery

os.environ.setdefault('DJANGO_SETTINGS_MODULE','celery_blast.settings')

app = Celery(
    'celery_blast')

#include=['blast_project.tasks','refseq_transactions.tasks','one_way_blast.tasks','external_tools.tasks']
app.config_from_object('django.conf:settings',namespace='CELERY')

app.autodiscover_tasks()

@app.task(bind=True)
def debug_task(self):
    print(f'Request: {self.request!r}')