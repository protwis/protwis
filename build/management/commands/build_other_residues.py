from django.core.management.base import BaseCommand, CommandError
from django.db.models import Q

from build.management.commands.build_human_residues import Command as BuildHumanResidues
from protein.models import ProteinConformation


class Command(BuildHumanResidues):
    help = 'Creates residue records for non-human receptors'

    pconfs = ProteinConformation.objects.filter(~Q(protein__species__id=1)).prefetch_related(
        'protein__residue_numbering_scheme__parent')
