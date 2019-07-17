import django
import os, sys

def setup_django():
    PWD = os.getenv("PWD")
    os.chdir(PWD)
    sys.path.insert(0, os.getenv("PWD"))
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "protwis.settings")
    django.setup()
