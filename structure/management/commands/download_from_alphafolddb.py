from django.core.management.base import BaseCommand
from django.conf import settings

from urllib.request import urlopen
import os


class Command(BaseCommand):
    help = '''Downloads AlphaFold models from EBI based on input UniProt accessions'''

    def add_arguments(self, parser):
        parser.add_argument('-u', '--uniprot',
            help='UniProt accession code(s)',
            nargs='+',
            default=False,
            type=str)
        parser.add_argument('-f' '--folder',
            dest='folder',
            action='store',
            type=str,
            default=False,
            help='Input folder of UniProt accession code named files')
        parser.add_argument('-o', '--output',
            dest='output',
            help='Folder to save files in')

    def handle(self, *args, **options):
        outf = options['output']
        if not os.path.exists(outf):
            os.mkdir(outf)

        accession_list = []

        if options['uniprot']:
            accession_list+=options['uniprot']
        if options['folder']:
            inputf = options['folder']
            files = os.listdir(inputf)
            for f in files:
                accession_list.append(f.split('.')[0])

        print(accession_list)
        for accession in accession_list:
            url = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'.format(accession)
            pdbdata_raw = urlopen(url).read().decode('utf-8')
            with open(os.sep.join([outf, accession+'.pdb']), 'w') as f:
                f.write(pdbdata_raw)
