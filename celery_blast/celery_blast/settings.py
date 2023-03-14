"""
    Django settings for celery_blast project.

    Generated by 'django-admin startproject' using Django 2.2.5.

    For more information on this file, see
    https://docs.djangoproject.com/en/2.2/topics/settings/

    For the full list of settings and their values, see
    https://docs.djangoproject.com/en/2.2/ref/settings/
"""

import os
from decouple import config

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/2.2/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = config('SECRET_KEY')

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = config("DEBUG", default=False, cast=bool)

ALLOWED_HOSTS = config('DJANGO_ALLOWED_HOSTS').split(" ")

# Application definition
INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_extensions',
    'blast_project',
    'refseq_transactions',
    'django_celery_results',
    'celery_progress',
    'one_way_blast',
    'external_tools',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'celery_blast.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': ['celery_blast/templates', 'static/images/result_images', 'media/blast_projects', 'media/one_way_blast',
                 'media/one_way_blast/remote_searches', 'media'],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'celery_blast.wsgi.application'

# Database
# https://docs.djangoproject.com/en/2.2/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': config('SQL_ENGINE', default="django.db.backends.sqlite3"),#'django.db.backends.postgresql_psycopg2',
        'NAME': config('SQL_DATABASE', default=BASE_DIR+'/db.sqlite3'),#'postgres',
        'USER': config('SQL_USER',default='user'),#'postgres',
        'PASSWORD': config('SQL_PASSWORD', default='password'),#'postgres',
        'HOST': config('SQL_HOST',default='localhost'),#'postgres',
        'PORT': config('SQL_PORT',default='5432')#'5432',
    }
}

# Password validation
# https://docs.djangoproject.com/en/2.2/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
# https://docs.djangoproject.com/en/2.2/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Berlin'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.2/howto/static-files/
STATIC_URL = '/static/'
STATIC_ROOT = '/blast/reciprocal_blast/assets'

MEDIA_URL = '/media/'
MEDIA_ROOT = '/media/'

STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "static"),
)

FILE_UPLOAD_TEMP_DIR = '/tmp'

# celery settings from celery projects website and following github repo: https://github.com/chrisk314/django-celery-docker-example/tree/master/mysite/settings
CELERY_BROKER_URL = 'pyamqp://rabbitmq:5672'
CELERY_RESULT_BACKEND = 'django-db'
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'

CELERY_TASK_SOFT_TIME_LIMIT = 345600  # 2 * 60min = 120min * 60sec = 7200sec
CELERY_TASK_TIME_LIMIT = 345660
SUBPROCESS_TIME_LIMIT = CELERY_TASK_SOFT_TIME_LIMIT - 5

'''
# django setting.
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'my_cache_table',
    }
}
'''

DATA_UPLOAD_MAX_MEMORY_SIZE = 5000000000
FILE_UPLOAD_MAX_MEMORY_SIZE = 5000000000

PANOPTES_IP = 'http://panoptes:5000'
STATIC_RESULT_IMAGES = STATIC_URL + 'images/result_images/'
BLAST_PROJECT_DIR = 'media/blast_projects/'
ONE_WAY_BLAST_PROJECT_DIR = 'media/one_way_blast/'
BLAST_DATABASE_DIR = 'media/databases/'
CDD_DIR = 'media/databases/CDD/Cdd'
ESEARCH_OUTPUT = 'media/esearch_output/'
CDD_DATABASE_URL = "https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz"
TAXDB_URL = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"
REFSEQ_ASSEMBLY_FILE = 'media/databases/refseq_summary_file/'
REFSEQ_URL = "ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
TAXONOMIC_NODES = 'media/taxonomic_node_files/'