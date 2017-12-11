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
    from protwis.settings_local_development import *

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
    'django.contrib.humanize',
    'debug_toolbar',
    'rest_framework',
    'rest_framework_swagger',
    'polymorphic',
    'common',
    'api',
    'news',
    'pages',
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
    'phylogenetic_trees',
    'sitesearch',
    'build_' + SITE_NAME,
    'construct',
    'tools',
    'drugs',
    'signprot',
    'mutational_landscape',
    'contactnetwork',
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
STATIC_ROOT = '/protwis/static/protwis'
STATICFILES_DIRS = (os.sep.join([BASE_DIR, "static"]),)
MEDIA_URL = '/media/'
MEDIA_ROOT = '/protwis/media/protwis'


# Serializer

SESSION_SERIALIZER = 'django.contrib.sessions.serializers.PickleSerializer'
SESSION_COOKIE_AGE = 86400 #Expire cookies and session after 24 hrs

SWAGGER_SETTINGS = {
    'USE_SESSION_AUTH' : False,
}

# Templates

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            # insert your TEMPLATE_DIRS here

        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.contrib.messages.context_processors.messages',
                'protwis.context_processors.google_analytics'
            ],
        },
    },
]

if DEBUG:
    TEMPLATES[0]['OPTIONS']['debug'] = True


# Debug toolbar
if DEBUG:
    DEBUG_TOOLBAR_PATCH_SETTINGS = False
    INTERNAL_IPS = ('10.0.2.2')


# Logging
if DEBUG:
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

#CACHE
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
        'LOCATION': '/tmp/django_cache',
        'OPTIONS': {
            'MAX_ENTRIES': 100000
        }
    },
    'alignments': {
        'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
        'LOCATION': '/tmp/django_cache_alignments',
        'OPTIONS': {
            'MAX_ENTRIES': 1000
        }
    }
}
