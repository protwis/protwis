# local settings
# override these settings by copying and editing this file to protwis/settings_local.py

# Site specific constants
SITE_NAME = 'gpcr' # used for site specific files
SITE_TITLE = 'GPCRdb' # for display in templates
ANALYTICS_KEY = ''
DATA_DIR = '/vagrant/data/protwis/' + SITE_NAME
DEFAULT_NUMBERING_SCHEME = 'gpcrdb'
DEFAULT_PROTEIN_STATE = 'inactive'
COMPARISON_SEGMENTS = ('TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7')


# Database

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'protwis',
        'USER': 'protwis',
        'PASSWORD': 'protwis',
        'HOST': 'localhost',
    }
}


# Quick-start development settings - unsuitable for production

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '-eqrx61@n*z3y1mc1w_@x1+yo(@^!k7i-vjaz0tx1$902a!4mu'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

# Debug toolbar
if DEBUG:
    DEBUG_TOOLBAR_PATCH_SETTINGS = False
    INTERNAL_IPS = ('10.0.2.2')

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []