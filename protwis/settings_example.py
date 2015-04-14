"""
Django settings for protwis project.

For more information on this file, see
https://docs.djangoproject.com/en/1.7/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.7/ref/settings/
"""

# Site specific variables
SITE_NAME = 'gpcr' # used for site specific files
SITE_TITLE = 'GPCRDB' # for display in templates
ANALYTICS_KEY = ''
DATA_DIR = '/vagrant/data/protwis/' + SITE_NAME
DEFAULT_NUMBERING_SCHEME = 'gpcrdb'


# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.7/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '-eqrx61@n*z3y1mc1w_@x1+yo(@^!k7i-vjaz0tx1$902a!4mu'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []


# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    # 'tastypie',
    'common',
    'home',
    'protein',
    'residue',
    'alignment',
    'similaritysearch',
    'build_' + SITE_NAME,
)

MIDDLEWARE_CLASSES = (
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


# Database
# https://docs.djangoproject.com/en/1.7/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'protwis',
        'USER': 'protwis',
        'PASSWORD': 'protwis',
        'HOST': 'localhost',
    }
}

# Internationalization
# https://docs.djangoproject.com/en/1.7/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Copenhagen'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.7/howto/static-files/

STATIC_URL = '/static/'
STATIC_ROOT = '/vagrant/protwis/static_root'
STATICFILES_DIRS = (os.sep.join([BASE_DIR, "static"]),)
MEDIA_URL = '/media/'
MEDIA_ROOT = '/vagrant/protwis/media'

# Serializer

SESSION_SERIALIZER = 'django.contrib.sessions.serializers.PickleSerializer'

# Logging

# LOGGING = {
#     'version': 1,
#     'disable_existing_loggers': False,
#     'formatters': {
#         'verbose': {
#             'format' : "[%(asctime)s] %(levelname)s [%(name)s:%(lineno)s] %(message)s",
#             'datefmt' : "%d/%b/%Y %H:%M:%S"
#         },
#     },
#     'handlers': {
#         'django': {
#             'level': 'DEBUG',
#             'class': 'logging.FileHandler',
#             'filename': 'logs/django.log',
#             'formatter': 'verbose'
#         },
#         'build': {
#             'level': 'DEBUG',
#             'class': 'logging.FileHandler',
#             'filename': 'logs/protwis_build.log',
#             'formatter': 'verbose'
#         },
#         'protwis': {
#             'level': 'DEBUG',
#             'class': 'logging.FileHandler',
#             'filename': 'logs/protwis.log',
#             'formatter': 'verbose'
#         },
#     },
#     'loggers': {
#         'django': {
#             'handlers':['django'],
#             'propagate': True,
#             'level':'DEBUG',
#         },
#         'build': {
#             'handlers': ['build'],
#             'level': 'DEBUG',
#         },
#         'protwis': {
#             'handlers': ['protwis'],
#             'level': 'DEBUG',
#         },
#     }
# }