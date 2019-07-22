from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from signprot.models import SignprotComplex

from contactnetwork.cube import *

import logging, json, os

class Command(BaseCommand):

    help = "Function to calculate interaction for only structures of GPCRs in complex with a signaling protein."

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):

        for pdb in SignprotComplex.objects.values_list('structure__pdb_code__index', flat=True):
            compute_interactions(pdb, True)
