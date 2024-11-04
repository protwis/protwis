# Local settings for django server running in Microsoft Windows
# WSL + Docker with PostgreSQL.
# Copy and edit the contents of this file to protwis/settings_local.py
import socket

# importing defaults
from protwis.settings_local_development import *
from protwis.settings_local_development import DATABASES

# Database. f-string requires Python 3.6+
DATABASES['default']['HOST'] = f'{socket.gethostname()}.local'


