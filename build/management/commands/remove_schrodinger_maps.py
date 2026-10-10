"""Remove the maps a build made for the Schrodinger imports, once both imports are done.

build_all builds chainmap.tsv (Engine 1 tree) and peptide_map.tsv (Engine 2
tree), imports with them, and then runs this, so the deliveries (a gpcrdb_data
checkout) are left as they came. It removes those two file names one level
below each tree, and a structure directory only when the map it just removed
was the last thing in it. Nothing else is touched, symlinked directories
included. When an import fails the caller stops before this step and the maps
stay for inspection.

    python manage.py remove_schrodinger_maps --engine1-dir <tree> --engine2-dir <tree>
"""

import os

from django.core.management.base import BaseCommand, CommandError

from interaction import schrodinger_import as si
from interaction import schrodinger_peptide as sp


def remove_maps(tree, name):
    """(maps removed, empty directories removed) under one tree."""
    maps = dirs = 0
    for entry in sorted(os.listdir(tree)):
        sub = os.path.join(tree, entry)
        if not os.path.isdir(sub) or os.path.islink(sub):
            continue
        path = os.path.join(sub, name)
        if not os.path.isfile(path):
            continue
        os.remove(path)
        maps += 1
        if not os.listdir(sub):
            os.rmdir(sub)
            dirs += 1
    return maps, dirs


class Command(BaseCommand):
    help = "Remove the chainmap.tsv and peptide_map.tsv files a build made, and the empty directories."

    def add_arguments(self, parser):
        parser.add_argument(
            "--engine1-dir", required=True, help="Engine 1 tree (chainmap.tsv)."
        )
        parser.add_argument(
            "--engine2-dir", required=True, help="Engine 2 tree (peptide_map.tsv)."
        )

    def handle(self, *args, **opt):
        for label, tree in (
            ("--engine1-dir", opt["engine1_dir"]),
            ("--engine2-dir", opt["engine2_dir"]),
        ):
            if not os.path.isdir(tree):
                raise CommandError("{} {!r} is not a directory".format(label, tree))
        for tree, name in (
            (opt["engine1_dir"], si.CHAINMAP_NAME),
            (opt["engine2_dir"], sp.MAP_NAME),
        ):
            maps, dirs = remove_maps(tree, name)
            self.stdout.write(
                "{}: removed {} {} and {} empty director{}".format(
                    tree, maps, name, dirs, "y" if dirs == 1 else "ies"
                )
            )
