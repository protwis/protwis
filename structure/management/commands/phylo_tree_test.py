from django.core.management.base import BaseCommand
from protein.models import ProteinFamily
from common import phylogenetic_tree as pt

from collections import OrderedDict
import json


class Command(BaseCommand):

    def handle(self, *args, **options):
        tree = pt.PhylogeneticTreeGenerator()
        qq = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
        #print('\n'.join([x[1].name for x in qq.tree.children.items()]))
        #print(qq.get_nodes(3))
        #print(json.dumps(qq.get_nodes_dict('crystals'), indent=4))
             
        #print(json.dumps(qq.get_nodes('crystalized')))
