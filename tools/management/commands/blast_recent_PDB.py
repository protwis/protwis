from Bio.Blast import NCBIXML
from django.core.management.base import BaseCommand
from django.conf import settings

from datetime import date
from dateutil.relativedelta import relativedelta
from io import StringIO
import logging, os
import requests
from subprocess import Popen, PIPE

from structure.models import Structure

class Command(BaseCommand):

    help = "Function to run against the local ."

    logger = logging.getLogger(__name__)
    rcsb_search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    rcsb_fasta_url = "https://www.rcsb.org/fasta"
    to_skip = ['8TII', '8TIL', '8TIN', '8TIO', '8VJ9', '8VTI', '9J7I', '9UOK', '9UT9', '9UTC', '9W1Z', '9OQ2', '9OQ3', '9OQ4', '9OQ5', '9OQ6'] # excluded from annotation, don't re-flag as new

    def add_arguments(self, parser):
        parser.add_argument('-m', '--months',
            dest='months',
            help='Number of months to check back',
            type=int,
            action='store',
            default=2)

    def handle(self, *args, **options):
        self.months = options['months']
        self.run()

    def run(self):
        # Request PDBs from RCSB for past X months
        match_date = (date.today() + relativedelta(months=-1*self.months)).strftime('%Y-%m-%d')

        rcsb_request = {"query":{
                            "type":"terminal",
                            "service":"text",
                            "parameters":{
                               "attribute":"rcsb_accession_info.initial_release_date",
                               "operator":"greater_or_equal",
                               "value": match_date + "T00:00:00Z"
                            }
                        },
                        "request_options":{
                            "return_all_hits": True
                        },
                       "return_type":"entry"
                    }
        rcsb_response = requests.post(self.rcsb_search_url, json=rcsb_request, headers={'Content-type': 'application/json'})

        if rcsb_response.status_code != 200:
            print("Incorrect response from RCSB web services - exiting")
            return

        # Collect all PDBs
        json_data = rcsb_response.json()
        new_pdbs = [entry["identifier"] for entry in json_data["result_set"]]
        print("Found", len(new_pdbs), "PDB entries")

        # Request sequences from RCSB (in sets of max pdbs per request)
        max_entries = 800 # it's actually 1000, but just to be on the safe side
        pdb_sets = [new_pdbs[i:i + max_entries] for i in range(0, len(new_pdbs), max_entries)]

        fasta_results = ""
        for pdb_set in pdb_sets:
            for header, sequence in grouped(self.fetch_fasta(pdb_set).splitlines(), 2):
                # Removal of RNA sequences and short sequences
                if "U" not in sequence and len(sequence) > 100:
                    fasta_results = fasta_results + header + "\n" + sequence + "\n"

        # BLAST against local BLAST database
        blast = Popen('%s -db %s -outfmt 5 -evalue 0.001 -max_target_seqs 1' % ('blastp',
            os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_human_bundle_blastdb'])),
            universal_newlines=True, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        (blast_out, blast_err) = blast.communicate(input=str(fasta_results))

        pdb_list = []
        if len(blast_err) != 0 and not blast_err.startswith('Warning'):
            print("BLAST search returned an error - exiting")
            return
        elif blast_out!='\n':
            # Process results and remove structures already present in GPCRdb
            blast_results = NCBIXML.parse(StringIO(blast_out))
            for result in blast_results:
                pdb_code = result.query.split('_')[0]
                if pdb_code in self.to_skip:
                    continue
                if len(result.alignments)>=1 and Structure.objects.filter(pdb_code__index=result.query[:4]).count() == 0:
                    top_hit = result.alignments[0].hsps[0]
                    if top_hit.score > 100:
                        print("HIT", "{0:>7}{1:>8}".format(top_hit.score, round(top_hit.expect,5)), result.query)
                        pdb_list.append(pdb_code)
        print(pdb_list)
        return pdb_list

    def fetch_fasta(self, pdb_ids):
        # RCSB's batch FASTA endpoint occasionally returns HTTP 200 with an error body
        # (e.g. "Error processing gql data: {}") for a single problematic entry in the
        # list, which would otherwise silently drop every entry in the batch. Bisect on
        # failure to isolate and skip only the offending id(s).
        post_data = {
                "structureIdList": ",".join(pdb_ids),
                "type": "entry",
                "outputType": "single"
            }
        rcsb_response = requests.post(self.rcsb_fasta_url, data=post_data)

        if rcsb_response.status_code == 200 and rcsb_response.text.strip().startswith(">"):
            return rcsb_response.text

        if len(pdb_ids) == 1:
            print("Skipping", pdb_ids[0], "- FASTA download failed")
            return ""

        mid = len(pdb_ids) // 2
        return self.fetch_fasta(pdb_ids[:mid]) + self.fetch_fasta(pdb_ids[mid:])

def grouped(iterable, n):
    return zip(*[iter(iterable)]*n)
