"""
Django settings for protwis project.
"""

# import local settings
# by default, local settings are in protwis/settings_local_default.py
# you can override these settings by creating a protwis/settings_local.py file (or copying settings_local_default)
# protwis/settings_local.py is ignored by git
try:
    from protwis.settings_local import *
except ImportError:
    from protwis.settings_local_default import *

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))


# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'debug_toolbar',
    'rest_framework',
    'rest_framework_swagger',
    'django_nvd3',
    'common',
    'api',
    'home',
    'protein',
    'family',
    'residue',
    'alignment',
    'similaritysearch',
    'similaritymatrix',
    'structure',
    'ligand',
    'interaction',
    'mutation',
    'build_' + SITE_NAME,
)

MIDDLEWARE_CLASSES = (
    'debug_toolbar.middleware.DebugToolbarMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'protwis.urls'

WSGI_APPLICATION = 'protwis.wsgi.application'


# Internationalization

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Copenhagen'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)

STATIC_URL = '/static/'
STATIC_ROOT = '/vagrant/protwis/static_root'
STATICFILES_DIRS = (os.sep.join([BASE_DIR, "static"]),)
MEDIA_URL = '/media/'
MEDIA_ROOT = '/vagrant/protwis/media'


# Serializer

SESSION_SERIALIZER = 'django.contrib.sessions.serializers.PickleSerializer'


# Logging

LOGGING = {
   'version': 1,
   'disable_existing_loggers': False,
   'formatters': {
       'verbose': {
           'format' : "[%(asctime)s] %(levelname)s [%(name)s:%(lineno)s] %(message)s",
           'datefmt' : "%d/%b/%Y %H:%M:%S"
       },
   },
   'handlers': {
       'django': {
           'level': 'DEBUG',
           'class': 'logging.FileHandler',
           'filename': 'logs/django.log',
           'formatter': 'verbose'
       },
       'build': {
           'level': 'DEBUG',
           'class': 'logging.FileHandler',
           'filename': 'logs/build.log',
           'formatter': 'verbose'
       },
       'protwis': {
           'level': 'DEBUG',
           'class': 'logging.FileHandler',
           'filename': 'logs/protwis.log',
           'formatter': 'verbose'
       },
   },
   'loggers': {
       'django': {
           'handlers':['django'],
           'propagate': True,
           'level':'DEBUG',
       },
       'build': {
           'handlers': ['build'],
           'level': 'DEBUG',
       },
       'protwis': {
           'handlers': ['protwis'],
           'level': 'DEBUG',
       },
   }
}