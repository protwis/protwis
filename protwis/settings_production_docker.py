# flake8: noqa: F401, F403, F405
"""
Docker production settings.

Use as: DJANGO_SETTINGS_MODULE=protwis.settings_production_docker

Inherits all site-specific constants (SITE_NAME, REFERENCE_POSITIONS,
HUB_ENABLED, etc.) from settings_local_production.py. This file only contains
what differs for the docker runtime: the database host, the data directory,
and the secrets/hosts that must come from the environment.

SECRET_KEY and DJANGO_ALLOWED_HOSTS are required and fail loudly if unset.
"""

import os
import sys

from protwis.settings_local_production import *  # noqa: F401, F403

DATA_DIR = os.environ.get("PROTWIS_DATA_DIR", "/protwis/data/protwis/" + SITE_NAME)
BUILD_CACHE_DIR = DATA_DIR + "/cache"

DATABASES["default"]["HOST"] = os.environ.get("POSTGRES_HOST", "db")
DATABASES["default"]["PORT"] = os.environ.get("POSTGRES_PORT", "5432")
DATABASES["default"]["NAME"] = os.environ.get(
    "POSTGRES_DB", DATABASES["default"]["NAME"]
)
DATABASES["default"]["USER"] = os.environ.get(
    "POSTGRES_USER", DATABASES["default"]["USER"]
)
DATABASES["default"]["PASSWORD"] = os.environ.get(
    "POSTGRES_PASSWORD", DATABASES["default"]["PASSWORD"]
)

# Required in production — fail loudly if either is missing.
SECRET_KEY = os.environ["DJANGO_SECRET_KEY"]
ALLOWED_HOSTS = [h for h in os.environ["DJANGO_ALLOWED_HOSTS"].split(",") if h]
DEBUG = False

# Register self as protwis.settings_local so settings.py picks up our values
# (notably DEBUG=False, so the if-DEBUG branches in settings.py do not fire).
sys.modules.setdefault("protwis.settings_local", sys.modules[__name__])

from protwis.settings import *  # noqa: E402, F403
