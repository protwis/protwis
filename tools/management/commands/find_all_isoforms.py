from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from protein.models import *
from common.tools import fetch_from_web_api
from Bio import pairwise2

import time
import json

from operator import itemgetter
from itertools import groupby
from collections import OrderedDict

class Command(BaseCommand):

    help = "Find all isoforms for GPCRs"

    def handle(self, *args, **options):

        ## Url API to map genename to ensemble ID
        cache_dir_genes = ['genenames', 'gene_lookup']
        url_gene = 'http://rest.genenames.org/fetch/symbol/$index'

        ensembl_version = 'grch37' # anything uses newest

        if ensembl_version=='grch37':
            ## Url to lookup ensemble ID to find transcripts
            cache_dir_transcripts = ['ensembl37', 'transcripts']
            url_ensembl = 'https://grch37.rest.ensembl.org/lookup/id/$index?expand=1;content-type=application/json'

            ## Url to lookup sequence of transcript
            cache_dir_seq = ['ensembl37', 'seq']
            url_ensembl_seq = 'https://grch37.rest.ensembl.org/sequence/id/$index?content-type=application/json'
        else:
            ## Url to lookup ensemble ID to find transcripts
            cache_dir_transcripts = ['ensembl', 'transcripts']
            url_ensembl = 'https://rest.ensembl.org/lookup/id/$index?expand=1;content-type=application/json'

            ## Url to lookup sequence of transcript
            cache_dir_seq = ['ensembl', 'seq']
            url_ensembl_seq = 'https://rest.ensembl.org/sequence/id/$index?content-type=application/json'
        

        # Get all human GPCRs
        ps = Protein.objects.filter(sequence_type__slug='wt', species__common_name="Human", family__slug__startswith='00').all().prefetch_related('genes').order_by('entry_name')
       
        isoforms = {}
        total_transcripts = 0
        total_proteins_with_isoforms = 0
        gene_to_ensembl = {}
        for p in ps:
            transcripts = []
            genes = list(p.genes.all().values_list('name',flat=True))
            print(">" + p.entry_name, 'genes:',genes)
            for gene in genes:

                # Use requests method due to weird functionality of genenames.org
                import requests
                url = 'http://rest.genenames.org/fetch/symbol/{}'.format(gene)
                cache_file_path = '{}/{}'.format('/'.join(cache_dir_genes), gene)
                # try fetching from cache
                data = cache.get(cache_file_path)
                if not data:
                    headers = {'Accept': 'application/json'}
                    try:
                        resp = requests.get(url=url, headers=headers)
                        data = resp.json()
                        cache.set(cache_file_path, data, 60*60*24*7) #7 days
                    except:
                        print('Error converting',gene)
                        continue
                if data['response']['docs']:
                    try:
                        # Get ensemble_gene_id
                        ensembl_gene_id = data['response']['docs'][0]['ensembl_gene_id']
                        gene_to_ensembl[p.entry_name] = ensembl_gene_id
                        #print("E_ID: " +ensembl_gene_id)
                        ensembl_transcripts = fetch_from_web_api(url_ensembl, ensembl_gene_id, cache_dir_transcripts)
                        for t in ensembl_transcripts['Transcript']:
                            display_name = t['display_name']
                            is_canonical = t['is_canonical']
                            if is_canonical:
                                # Skip canonical entries
                                continue
                            biotype = t['biotype']
                            t_id = t['id']

                            # Only interested in protein_coding
                            if biotype=='protein_coding':
                                length = t['Translation']['length']
                                seq_id = t['Translation']['id']
                                transcript_info = OrderedDict([('display_name',display_name),('t_id',t_id),('length',length), ('seq_id',seq_id)])
                                seq = fetch_from_web_api(url_ensembl_seq, seq_id,cache_dir_seq)
                                transcript_info['seq'] = seq['seq']
                                transcripts.append(transcript_info)
                                total_transcripts += 1
                    except:
                        print('Error fetching ensemble_gene_id for gene',gene)
                        pass
            print(len(transcripts), 'transcripts found')
            
            # Add if transcripts found
            if len(transcripts):
                isoforms[p.entry_name] = transcripts
                total_proteins_with_isoforms += 1

        # print small summary results
        print('total_proteins_searched',len(ps))
        print('total_proteins_with_isoforms', total_proteins_with_isoforms)
        print('total_transcripts',total_transcripts)

        print(gene_to_ensembl)
        # save to file
        f = open('protein/data/all_isoforms.json', 'w')
        json.dump(isoforms,f, indent=4, separators=(',', ': '))