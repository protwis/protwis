from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection

from structure.models import Structure

import logging, json, os
import time

class Command(BaseCommand):

    help = "Function evaluating sequence compatibility with wildtype protein."

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        Structure.objects.filter(refined=True).delete()
