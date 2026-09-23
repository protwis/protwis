# flake8: noqa: F401, F403, F405
"""
Docker development settings.

Use as: DJANGO_SETTINGS_MODULE=protwis.settings_local_docker

Inherits all site-specific constants (SITE_NAME, REFERENCE_POSITIONS,
DOCUMENTATION_URL, dummy CACHES, etc.) from settings_local_development.py.
This file only contains what differs for the docker runtime: the database
host, the data directory, and any value supplied via environment variables.

Single source of truth: edit shared constants in settings_local_development.py
and this file picks up the change automatically.
"""

import os
import sys

from protwis.settings_local_development import *  # noqa: F401, F403

# GPCRdb data path inside the container (bind-mounted by docker-compose).
DATA_DIR = os.environ.get("PROTWIS_DATA_DIR", "/app/data/protwis/" + SITE_NAME)
BUILD_CACHE_DIR = DATA_DIR + "/cache"

# Point at the docker-compose `db` service by default; everything else
# falls back to whatever settings_local_development.py defined.
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

# Optional env overrides for SECRET_KEY / DEBUG / ALLOWED_HOSTS.
SECRET_KEY = os.environ.get("DJANGO_SECRET_KEY", SECRET_KEY)
DEBUG = os.environ.get("DJANGO_DEBUG", str(DEBUG)).lower() in ("1", "true", "yes", "on")
ALLOWED_HOSTS = os.environ.get("DJANGO_ALLOWED_HOSTS", ",".join(ALLOWED_HOSTS)).split(
    ","
)

# Register self as protwis.settings_local so settings.py's
# `try: from protwis.settings_local import *` picks up our env-driven values.
# Lets us reuse all of base settings.py without duplicating it.
sys.modules.setdefault("protwis.settings_local", sys.modules[__name__])

from protwis.settings import *  # noqa: E402, F403

# Re-apply dummy caches: settings.py defines a FileBasedCache CACHES near the
# bottom that overwrites the dummy CACHES inherited from settings_local_development.
CACHES = {
    "default": {"BACKEND": "django.core.cache.backends.dummy.DummyCache"},
    "alignments": {"BACKEND": "django.core.cache.backends.dummy.DummyCache"},
}
