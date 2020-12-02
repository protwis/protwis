from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from protein.models import *
from common.tools import fetch_from_web_api
from Bio import pairwise2

from structure.functions import BlastSearch

import time
import json
import os

from operator import itemgetter
from itertools import groupby
from collections import OrderedDict, defaultdict

class Command(BaseCommand):

    help = "Find all isoform sequences"

    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'protein_data', 'uniprot'])
    remote_uniprot_dir = 'http://www.uniprot.org/uniprot/'

    def parse_uniprot_file(self, accession):

        filename = accession + '.txt'
        local_file_path = os.sep.join([self.local_uniprot_dir, filename])
        remote_file_path = self.remote_uniprot_dir + filename

        up = {}
        up['genes'] = set()
        up['proteins'] = set()
        up['transcripts'] = set()

        read_sequence = False
        remote = False

        # record whether organism has been read
        os_read = False

        # should local file be written?
        local_file = False

        #try:
        if os.path.isfile(local_file_path):
            uf = open(local_file_path, 'r')
        else:
            uf = urlopen(remote_file_path)
            remote = True
            local_file = open(local_file_path, 'w')

        for raw_line in uf:
            # line format
            if remote:
                line = raw_line.decode('UTF-8')
            else:
                line = raw_line

            # write to local file if appropriate
            if local_file:
                local_file.write(line)

            # end of file
            if line.startswith('//'):
                break

            # entry name and review status
            if line.startswith('DR'):
                split = line.split(";")
                if line.split()[1]=='Ensembl;':
                    gene = split[-1].split(".")[0].strip()
                    protein = split[2].strip()
                    transcript = split[1].strip()
                    up['genes'].add(gene)
                    up['proteins'].add(protein)
                    up['transcripts'].add(transcript)

        # close the Uniprot file
        uf.close()
        # except:
        #     print('uniprot failed!')
        #     return False

        # close the local file if appropriate
        if local_file:
            local_file.close()

        return up

    def find_ensembl_id(self,gene):
        # Use requests method due to weird functionality of genenames.org
        import requests
        cache_dir_genes = ['genenames', 'gene_lookup']
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
                return False
        if data['response']['docs']:
            try:
                # Get ensemble_gene_id
                return data['response']['docs'][0]['ensembl_gene_id']
            except:
                return False
        else:
            return False

    def find_ensembl_id_by_uniprot(self,uniprot_id):

        uniprot = self.parse_uniprot_file(uniprot_id)
        for k,v in uniprot.items():
            uniprot[k] = list(v)
        return uniprot

    def handle(self, *args, **options):

        ## Prepare comparasion info ##
        filepath = 'protein/data/Isoform_annotation_table.txt'
        lmb_data = OrderedDict()
        total_lmb_isoforms = 0
        all_lmb_isoforms = []
        with open(filepath, "r", encoding='UTF-8') as f:
            for i,row in enumerate(f):
                if i>0:
                    c = row.split("\t")
                    entry_name = "{}_human".format(c[1].lower())
                    transcripts = c[4].split(", ")

                    if not entry_name in lmb_data:
                        lmb_data[entry_name] = []
                    lmb_data[entry_name] += transcripts
                    total_lmb_isoforms += 1
                    all_lmb_isoforms += transcripts

        print('all_lmb_isoforms',len(all_lmb_isoforms),'distinct',len(set(all_lmb_isoforms)))

        
        ## Get parsed gtex annotation
        with open('protein/data/matched_gtex.json') as json_file:
            gtex_old = json.load(json_file)

        ## Need to rewrite these entries, as ensembl doesnt use the . for transcripts
        gtex = {}
        for key, val in gtex_old['transcripts'].items():
            t,g = key.split("_")
            new_key = "{}_{}".format(t.split(".")[0],g)
            gtex[new_key] = val
            # del gtex[new_key]['subjects']


        ## Url API to map genename to ensemble ID
        cache_dir_genes = ['gtexportal', 'gene_lookup']
        url_gene = 'https://gtexportal.org/rest/v1/reference/gene?geneId=$index&gencodeVersion=v19&genomeBuild=GRCh37%2Fhg19&pageSize=250&format=json'

        ## Url to lookup ensemble ID to find transcripts
        cache_dir_transcripts_gtex = ['gtexportal', 'transcripts']
        url_transcripts = 'https://gtexportal.org/rest/v1/reference/transcript?gencodeId=$index&gencodeVersion=v19&genomeBuild=GRCh37%2Fhg19'

        cache_dir_transcripts = ['ensembl37', 'transcripts']
        url_ensembl = 'https://grch37.rest.ensembl.org/lookup/id/$index?expand=1;content-type=application/json'

        cache_dir_gtex_expression  = ['gtexportal', 'expression_data']
        url_expression = 'https://gtexportal.org/rest/v1/expression/medianTranscriptExpression?datasetId=gtex_v7&gencodeId=$index&format=json'

        ## Url to lookup sequence of transcript
        cache_dir_seq = ['ensembl37', 'seq_protein']
        url_ensembl_seq = 'https://grch37.rest.ensembl.org/sequence/id/$index?content-type=application/json;type=protein'

        # Get all human GPCRs
        ps = Protein.objects.filter(sequence_type__slug='wt', species__common_name="Human", family__slug__startswith='00').all().prefetch_related('genes').order_by('entry_name')
       
        isoforms = {}
        total_transcripts = 0
        total_transcript_skipped_no_tissue=0
        total_proteins_with_isoforms = 0
        gene_to_ensembl = {}
        transcripts_ids_total = set()
        transcripts_ids_skipped_total = set()
        total_fetched_transcripts = 0
        canonical_disagreement_count = 0

        total_new_transcripts = []
        total_not_found = []
        total_not_found_due_to_skipped = []
        new_proteins = set()

        lmb_compare_sequences = [0,0,0] # correct, wrong, not exists in lmb

        sequence_lookup = {}

        ## COMPARE SEQUENCES
        filenames = os.listdir("protein/data/LMB_sequences/")
        all_lmb_sequences= {}
        for f in filenames:
            with open ("protein/data/LMB_sequences/"+f, "r") as myfile:
                fasta=myfile.read().splitlines()
                for i,l in enumerate(fasta):
                    if l[0]==">":
                        e_id = l[2:]
                        continue
                    if e_id in all_lmb_sequences:
                        print('already there!',e_id)
                    if i>2:
                        all_lmb_sequences[e_id]=l
        print('all_lmb_sequences',len(all_lmb_sequences))

        f = open("protein/data/20190726_transcripts.fa", "w")
        missing_sequences = 0
        total_lmb_sequences = 0
        sequences_lookup = defaultdict(list)
        for p,ts in lmb_data.items():
            seq = Protein.objects.get(entry_name=p).sequence
            sequences_lookup[seq].append([p,p])
            # print(p,ts)
            # print(seq)
            f.write(">{} GPCRdb sequence reference\n".format(p))
            f.write("{}\n".format(seq))
            seq_filename = "protein/data/LMB_sequences/{}_nonstrict_transcripts.fa".format(p)
            lmb_sequences = {}
            try:
                with open (seq_filename, "r") as myfile:
                    #fasta_raw = myfile.read()
                    fasta=myfile.read().splitlines()
                    for i,l in enumerate(fasta):
                        if l[0]==">":
                            e_id = l[2:]
                            continue
                        lmb_sequences[e_id]=l
                        if i>2:
                            total_lmb_sequences += 1
            except:
                #print('No file for',p,' So no sequence for',ts)
                missing_sequences += len(ts)
            for t in ts:
                if not t in lmb_sequences:
                    #print('missing ',t,'in',"{}_nonstrict_transcripts.fa".format(p))
                    missing_sequences += 1

                seq = fetch_from_web_api(url_ensembl_seq, t,cache_dir_seq)['seq']
                sequences_lookup[seq].append([t,p])
                if t in lmb_sequences:
                    if seq!=lmb_sequences[t]:
                        print(t,'different from LBM - length ensembl:',len(seq),"length lmb:",len(lmb_sequences[t]))
                f.write(">{} ({})\n".format(t,p))
                f.write("{}\n".format(seq))
        f.close()
        print('total missing sequences',missing_sequences)
        print('total lmb transcript sequences provided',total_lmb_sequences)
        print('total lmb protein',len(lmb_data))
        #return
        for seq,ts in sequences_lookup.items():
            if len(ts)>1:
                print('Identical sequence:',ts)

        sequences_lookup = defaultdict(list)
        all_transcript_seq = {}
        for p in ps:# .filter(entry_name='gpc5b_human').all():
            transcripts = []
            transcripts_ids = []
            transcripts_ids_skipped = []
            ensembl_transcripts_count = 0
            genes = list(p.genes.all().values_list('name',flat=True))
            uniprot = p.accession
            canonical = ''
            canon_seq = p.sequence
            # sequence_lookup[canon_seq] = p.entry_name
            grch37_canonical_seq = ''
            uniprot_canonical = ''
            grch37_canonical = ''

            # print(">" + p.entry_name,uniprot, 'genes:',genes)
            seq_filename = "protein/data/LMB_sequences/{}_nonstrict_transcripts.fa".format(p.entry_name)
            lmb_sequences = {}
            try:
                with open (seq_filename, "r") as myfile:
                    #fasta_raw = myfile.read()
                    fasta=myfile.read().splitlines()
                    for l in fasta:
                        if l[0]==">":
                            e_id = l[2:]
                            continue
                        lmb_sequences[e_id]=l
            except:
                pass
            #break

            alternative_ids_uniprot = self.find_ensembl_id_by_uniprot(uniprot)
            # print(alternative_ids_uniprot)
            ensembl_gene_id = []
            for gene in genes:
                if not gene:
                    continue
                gene_lookup = fetch_from_web_api(url_gene, gene, cache_dir_genes)
                
                # try:
                same_gene_id = ''
                if gene_lookup['gene']:
                    for gene_info in gene_lookup['gene']:
                        if gene_info['geneSymbol']==gene:
                            ensembl_gene_id.append(gene_info['gencodeId'])

            if len(ensembl_gene_id)>1:
                print(ensembl_gene_id,'MORE THAN 1 !!!!')

            if len(ensembl_gene_id)==0:
                print('No ID found, using uniprot')
                if alternative_ids_uniprot['genes']:
                    ensembl_gene_id = alternative_ids_uniprot['genes'][0]
                else:
                    print("NO ID FOR THIS RECEPTOR")
                    continue
            else:
                ensembl_gene_id = ensembl_gene_id[0]

            #alternative_id = self.find_ensembl_id(gene)
            # alternative_id_uniprot = self.find_ensembl_id_by_uniprot(uniprot)
            # print(ensembl_gene_id,alternative_ids_uniprot)
            # expression = fetch_from_web_api(url_expression,ensembl_gene_id,cache_dir_gtex_expression)
            # print(expression)
            # go through expression
            # expressed_transcripts = {}
            # for e in expression['medianTranscriptExpression']:
            #     if e['median']>0 or 1==1:
            #         #only if expression
            #         t_id = e['transcriptId']
            #         t_short = t_id.split(".")[0]
            #         tissue = e['tissueSiteDetailId']
            #         if t_short not in expressed_transcripts:
            #             expressed_transcripts[t_short] = {'long':t_id,'tissues':[], 'max_median':0}
            #         if expressed_transcripts[t_short]['max_median']<e['median']:
            #             expressed_transcripts[t_short]['max_median'] = e['median']
            #         expressed_transcripts[t_short]['tissues'].append([tissue,e['median']])   
            # print(expressed_transcripts)
            # print(ensembl_gene_id)
            gene_to_ensembl[p.entry_name] = ensembl_gene_id
            # print("E_ID: " +ensembl_gene_id,alternative_ids_uniprot)
            # ensembl_transcripts = fetch_from_web_api(url_ensembl, ensembl_gene_id, cache_dir_transcripts)
            # use uniprot gene ID instead
            ensembl_transcripts = fetch_from_web_api(url_ensembl, ensembl_gene_id, cache_dir_transcripts)
            # print(ensembl_gene_id)
            if (alternative_ids_uniprot['genes'] and ensembl_gene_id.split(".")[0]!=alternative_ids_uniprot['genes'][0]):
                print("##### ensembl gene id changed",ensembl_gene_id,alternative_ids_uniprot['genes'][0])

            #total_fetched_transcripts += len(ensembl_transcripts['Transcript'])
            # print(ensembl_transcripts)
            same_gene_id = True
            if not ensembl_transcripts:
                print('error',alternative_ids_uniprot,ensembl_gene_id)
                same_gene_id = False
                ensembl_transcripts = fetch_from_web_api(url_ensembl, alternative_ids_uniprot['genes'][0], cache_dir_transcripts)

            for t in ensembl_transcripts['Transcript']:
                ensembl_transcripts_count += 1
                display_name = t['display_name']
                is_canonical = t['is_canonical']
                biotype = t['biotype']
                t_id = t['id']
                #     # Skip canonical entries
                #     continue

                # Only interested in protein_coding
                if biotype=='protein_coding':
                    total_fetched_transcripts += 1

                    key = '{}_{}'.format(t_id,ensembl_gene_id)

                    if not key in gtex:
                        # print('t_id', t_id, 'not in expressed_transcripts')
                        total_transcript_skipped_no_tissue += 1
                        transcripts_ids_skipped_total.add(t_id)
                        transcripts_ids_skipped.append(t_id)
                        continue

                    if gtex[key]["count"]<3:
                        total_transcript_skipped_no_tissue += 1
                        transcripts_ids_skipped_total.add(t_id)
                        transcripts_ids_skipped.append(t_id)
                        continue

                    length = t['Translation']['length']
                    seq_id = t['Translation']['id']
                    transcript_info = OrderedDict([('display_name',display_name),('t_id',t_id),('length',length), ('seq_id',seq_id), ('expressed',gtex[key])])
                    seq = fetch_from_web_api(url_ensembl_seq, seq_id,cache_dir_seq)


                    if is_canonical:
                        grch37_canonical = t_id
                        transcript_info['grch37_canonical'] = True
                        grch37_canonical_seq = seq['seq']

                    if seq['seq']==canon_seq:
                        uniprot_canonical = t_id
                        transcript_info['uniprot_canonical'] = True
                        continue
                        # Skip canonical entries

                    sequences_lookup[seq['seq']].append([t_id,p.entry_name])
                    all_transcript_seq[t_id] = seq['seq']
                    if seq['seq'] in sequence_lookup:
                        print('SEQUENCE ALREADY SEEN',t_id, sequence_lookup[seq['seq']])
                        continue
                    sequence_lookup[seq['seq']] = t_id


                    transcript_info['seq'] = seq['seq']
                    if not t_id in lmb_sequences:
                        transcript_info['lmb_sequences'] = False
                        lmb_compare_sequences[2] += 1
                    else:
                        if lmb_sequences[t_id]==seq['seq']:
                            transcript_info['lmb_sequences'] = True
                            lmb_compare_sequences[0] += 1
                        else:
                            transcript_info['lmb_sequences'] = lmb_sequences[t_id]
                            lmb_compare_sequences[1] += 1

                    if t_id in alternative_ids_uniprot['transcripts']:
                        transcript_info['in_uniprot'] = True
                    else:
                        transcript_info['in_uniprot'] = False

                    if p.entry_name in lmb_data and t_id in lmb_data[p.entry_name]:
                        transcript_info['in_lmb'] = True
                    else:
                        transcript_info['in_lmb'] = False

                    if t_id not in transcripts_ids:
                        transcripts.append(transcript_info)
                        transcripts_ids.append(t_id)
                        transcripts_ids_total.add(t_id)
                    total_transcripts += 1
                # except:
                #     print('Error fetching ensemble_gene_id for gene',gene)
                #     pass

            not_found = []
            not_found_due_to_skipped = []
            if p.entry_name in lmb_data:
                for t in lmb_data[p.entry_name]:
                    if t not in transcripts_ids:
                        if t in transcripts_ids_skipped:

                            f = open("protein/data/20190726_skipped_due_to_gtex.txt", "a")
                            not_found_due_to_skipped.append(t)
                            key = '{}_{}'.format(t,ensembl_gene_id)
                            if not key in gtex:
                                reason = 'Not in GTEX'
                            else:
                                reason = 'Subjects low in GTEX - count is {} - subject ids {}'.format(gtex[key]['count'],", ".join(gtex[key]['subjects']))
                            f.write("{}: {}\n".format(t,reason))
                            f.close()
                            # print(t)
                        else:
                            not_found.append(t)

            total_not_found += not_found
            total_not_found_due_to_skipped += not_found_due_to_skipped

            new = []
            for t in transcripts_ids:
                if p.entry_name in lmb_data and t in lmb_data[p.entry_name]:
                    pass
                else:
                    ts_check = sequences_lookup[all_transcript_seq[t]]
                    for t_check in ts_check:
                        if p.entry_name in lmb_data and t_check in lmb_data[p.entry_name]:
                            print('found via duplicate',t_check,t)
                            continue
                    key = '{}_{}'.format(t,ensembl_gene_id)



                    #blast = BlastSearch(top_results=2)

                    blast = BlastSearch(blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_human_blastdb']), top_results=2)
                    blast_out = blast.run(all_transcript_seq[t])
                    result = [(Protein.objects.get(pk=x[0]).entry_name, x[1].hsps[0].expect) for x in blast_out]
                    #print(result)
                    if result:
                        if result[0][0]==p.entry_name and result[0][1]<0.05:
                            f = open("protein/data/20190726_new_transcripts_for_consideration.txt", "a")
                            reason = 'GTEX count: {}'.format(gtex[key]['count'])
                            f.write(">{} ({}): {}\n".format(t,p.entry_name,reason))
                            f.write("{}\n".format(all_transcript_seq[t]))
                            f.close()
                            new.append(t)
                            if p.entry_name in lmb_data:
                                new_proteins.add(p.entry_name)
                        else:
                            print('bad blast match',result)
                    else:
                        print('bad blast match',result)

            total_new_transcripts += new

            # print(len(alternative_ids_uniprot['transcripts']), 'uniprot transcripts found',ensembl_transcripts_count, ' ensembl transcripts found',len(transcripts), 'transcripts kept after filtering')
            
            # Add if transcripts found
            if len(transcripts):
                isoforms[p.entry_name] = {'ensembl_gene_id':ensembl_gene_id,'same_gene_id':same_gene_id,'canonical_seq':canon_seq, 'grch37_canonical_seq':grch37_canonical_seq, 'isoforms': transcripts, 'uniprot_lookup': alternative_ids_uniprot, 'lmb_not_found':not_found, 'lmb_not_found_due_to_skipped': not_found_due_to_skipped, 'new_transcripts_than_lmb': new,'skipped_due_to_gtex': transcripts_ids_skipped, 'grch37_canonical':grch37_canonical, 'uniprot_canonical':uniprot_canonical}
                if grch37_canonical_seq!=canon_seq:
                    isoforms[p.entry_name]['canonical_disagreement'] = True
                    canonical_disagreement_count += 1
                # isoforms[p.entry_name].append(alternative_ids_uniprot)
                # isoforms[p.entry_name].append(not_found)
                total_proteins_with_isoforms += 1
            else:
                isoforms[p.entry_name] = {'ensembl_gene_id':ensembl_gene_id,'same_gene_id':same_gene_id,'canonical_seq':canon_seq, 'grch37_canonical_seq':grch37_canonical_seq, 'isoforms': transcripts, 'uniprot_lookup': alternative_ids_uniprot, 'lmb_not_found':not_found, 'lmb_not_found_due_to_skipped': not_found_due_to_skipped, 'new_transcripts_than_lmb': new,'skipped_due_to_gtex': transcripts_ids_skipped, 'grch37_canonical':grch37_canonical, 'uniprot_canonical':uniprot_canonical}
                
            # break
            f = open('protein/data/all_isoforms_gtex.json', 'w')
            json.dump(isoforms,f, indent=4, separators=(',', ': '))
            #break

        for seq,ts in sequences_lookup.items():
            if len(ts)>1:
                print('identical sequence',ts)

        for t in total_not_found:
            ts_check = sequences_lookup[all_transcript_seq[t]]
            found = False
            for t_check in ts_check:
                if t_check[0] not in total_not_found:
                    print(t,'found but under another id',t_check[0])
                    found = True

            if not found:
                print('##',t,'in LMB but not in this search')


        # print small summary results
        print('total_proteins_searched',len(ps))
        print('total_proteins_with_isoforms', total_proteins_with_isoforms)
        print('Total transcripts deemed to be isoforms',total_transcripts)
        print('Amount of these not in LMB data',len(total_new_transcripts))
        print(new_proteins)
        # print('Amount in LBM not found',len(total_not_found))
        # print(total_not_found)
        print('Amount in LBM found but skipped due to GTEX data',len(total_not_found_due_to_skipped))
        print(total_not_found_due_to_skipped)
        print('Sequence compare to LMB', lmb_compare_sequences)
        print('canonical_disagreement_count',canonical_disagreement_count)
        print(total_not_found)
        # print('total_transcript_skipped_no_tissue',total_transcript_skipped_no_tissue)
        # print('total_transcript_skipped_no_tissue2 ',len(transcripts_ids_skipped_total))
        # print('total_fetched_transcripts',total_fetched_transcripts)

        # print(gene_to_ensembl)
        # save to file
        f = open('protein/data/all_isoforms_gtex.json', 'w')
        json.dump(isoforms,f, indent=4, separators=(',', ': '))