# local settings
# override these settings by copying and editing this file to protwis/settings_local.py

# Site specific constants
SITE_NAME = 'gpcr' # used for site specific files
SITE_TITLE = 'GPCRdb' # for display in templates
DATA_DIR = '/vagrant/data/protwis/' + SITE_NAME
DEFAULT_NUMBERING_SCHEME = 'gpcrdb'
DEFAULT_PROTEIN_STATE = 'inactive'
REFERENCE_POSITIONS = {'TM1': '1x50', 'TM2': '2x50', 'TM3': '3x50', 'TM4': '4x50', 'TM5': '5x50', 'TM6': '6x50',
    'TM7': '7x50', 'H8': '8x50'}

# analytics
GOOGLE_ANALYTICS_KEY = False

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