from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from phylogenetic_trees.views import *

class Command(BaseCommand):
    classes=['001','002','003','004','006','007']

    def handle(self, *args, **options):
        for cl in self.classes:
            Tree = Treeclass()
            Tree.Prepare_file('',cl)


if __name__ == '__main__':
    z = Command()
    z.handle()
