# Local settings for django server running in Microsoft Windows 
# WSL + Docker with PostgreSQL.
# Copy and edit the contents of this file to protwis/settings_local.py
import socket

# importing defaults
from protwis.settings_local_development import *

# Database. f-string requires Python 3.6+
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'protwis',
        'USER': 'protwis',
        'PASSWORD': 'protwis',
        'HOST': f'{socket.gethostname()}.local',
    }
}

