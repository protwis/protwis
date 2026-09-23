from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

import os
import yaml


ORDERED_SEGMENTS = ['1', 'i1', '2', 'e1', '3', 'i2', '4', 'e2', '5', '6', '7', '8']


def to_num(v):
    if v is None or v == '-':
        return None
    return float(v)


def check_entry(name, seg):
    """Check that, walking ORDERED_SEGMENTS in order, each segment's begin <=
    x50-middle <= end, and each segment's begin is higher than the last
    non-missing segment end seen so far."""
    errors = []
    last_end = None
    last_end_label = None
    for label in ORDERED_SEGMENTS:
        b = to_num(seg.get(label + 'b', '-'))
        x = to_num(seg.get(label + 'x', '-'))
        e = to_num(seg.get(label + 'e', '-'))

        if x is not None and not (b is not None and e is not None and b <= x <= e):
            errors.append((name, label, 'within-segment ascension', b, x, e))

        if b is not None and last_end is not None and not (last_end < b):
            errors.append((name, label, 'not higher than {} end'.format(last_end_label), last_end, None, b))

        if e is not None:
            last_end, last_end_label = e, label

    return errors


class Command(BaseCommand):
    help = '''Checks non_xtal_segends.yaml or mod_xtal_segends.yaml for segend consistency: within each
              segment begin <= x50-middle <= end, and each segment's begin is higher than the end of the
              last non-missing segment seen so far, walking segments in order {}.'''.format(ORDERED_SEGMENTS)

    file_paths = {
        'non_xtal': os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'non_xtal_segends.yaml']),
        'mod_xtal': os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'mod_xtal_segends.yaml']),
    }

    def add_arguments(self, parser):
        parser.add_argument('--file',
            choices=['non_xtal', 'mod_xtal'],
            default='non_xtal',
            help='Which annotation file to check (default: non_xtal)')
        parser.add_argument('--path',
            default=None,
            help='Custom path to a yaml file with the same schema, overrides --file')
        parser.add_argument('--entry',
            nargs='+',
            default=None,
            help='Only check these top-level entries (protein entry_names or PDB codes)')
        parser.add_argument('--verbose',
            action='store_true',
            default=False,
            help='Print the actual b/x/e values alongside each violation')

    def handle(self, *args, **options):
        path = options['path'] or self.file_paths[options['file']]

        with open(path, 'r') as f:
            data = yaml.safe_load(f)

        entries = options['entry'] if options['entry'] else list(data.keys())

        all_errors = []
        checked = 0
        for name in entries:
            if name not in data:
                self.stderr.write('WARNING: entry {} not found in {}'.format(name, path))
                continue
            checked += 1
            all_errors.extend(check_entry(name, data[name]))

        for name, label, reason, b_or_last_end, x, e_or_b in all_errors:
            if x is not None:
                msg = 'ERROR: {} segment {} {}: {} <= {} <= {}'.format(name, label, reason, b_or_last_end, x, e_or_b)
            else:
                msg = 'ERROR: {} segment {} {}: {} < {}'.format(name, label, reason, b_or_last_end, e_or_b)
            if options['verbose']:
                self.stdout.write(msg)
            else:
                self.stdout.write('ERROR: {} segment {} {}'.format(name, label, reason))

        self.stdout.write('Checked {} entries in {}'.format(checked, path))
        self.stdout.write('{} violations found'.format(len(all_errors)))

        if all_errors:
            raise CommandError('{} segend consistency violations found in {}'.format(len(all_errors), path))
