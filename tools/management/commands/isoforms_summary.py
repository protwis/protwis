from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q
from django.template.loader import render_to_string
from protein.models import *
from common.tools import fetch_from_web_api
from Bio import pairwise2

import time
import json

from operator import itemgetter
from itertools import groupby

class Command(BaseCommand):

    help = "Test isoforms"


    def handle(self, *args, **options):
        cache_dir = ['ensembl2', 'isoform']
        url = 'https://rest.ensembl.org/sequence/id/$index?content-type=application/json&type=protein'
        #Altenrative
        url = 'https://grch37.rest.ensembl.org/sequence/id/$index?db_type=core;object_type=predictiontranscript;content-type=application/json;species=homo_sapiens;type=protein'
        url = 'https://grch37.rest.ensembl.org/sequence/id/$index?type=protein;content-type=application/json'
        filepath = 'protein/data/Isoform_annotation_table.txt'
        isoforms = []
        with open(filepath, "r", encoding='UTF-8') as f:
            for row in f:
                c = row.split("\t")
                isoforms.append(c)

        # Skip header
        total_matches = 0
        total_mismatches = 0
        total_mismatches_1 = 0
        total_align_match = 0
        total_align_mismatch = 0

        isoforms_with_issue = {}

        dump = {}

        for c, i in enumerate(isoforms[1:]):
            
            p = '{}_human'.format(i[0].lower())
            print(p)
            protein = Protein.objects.get(entry_name=p, sequence_type__slug='wt', species__common_name='Human')
            wt_seq = protein.sequence
            rs = Residue.objects.filter(protein_conformation__protein=protein).prefetch_related('protein_segment','display_generic_number','generic_number')
            r_lookup = {}
            r_segment = {}
            for r in rs:
                r_lookup[r.sequence_number] = [r.protein_segment.slug,str(r.display_generic_number), r.sequence_number]
                if r.protein_segment.slug not in r_segment:
                    r_segment[r.protein_segment.slug] = 0
                r_segment[r.protein_segment.slug] += 1

            seq_filename = "protein/data/MSA_GPCR_isoforms/{}_isoform_MSA.fa".format(p.lower())
            with open (seq_filename, "r") as myfile:
                fasta_raw = myfile.read()
                fasta=fasta_raw.splitlines() 


            wt_seq2=fasta[1]

            es = i[3].split(", ")
            isoform_id = i[1]

            print(c, len(isoforms), p,isoform_id,es)

            wt_check = wt_seq==wt_seq2.replace("-","")
            if not wt_check:
                print(p, 'WT SEQ NO MATCH!!')
                # continue
            # print('WT SEQ',wt_seq==wt_seq2.replace("-",""))
            ranges = {}
            for e in es:
                iso_seq_msa = fasta[1+int(isoform_id)*2]
                iso_seq_msa_corrected = ''
                for pos, a in enumerate(iso_seq_msa):
                    if wt_seq2[pos]=='-' and a=='-':
                        continue
                    iso_seq_msa_corrected += a

                isoform_info = fetch_from_web_api(url, e, cache_dir)
                if (isoform_info):
                    iso_seq = isoform_info['seq']
                    iso_check = iso_seq == iso_seq_msa.replace("-","")
                    if not iso_check:
                        isoforms_with_issue[p+"_"+e] = "Sequence does not match with API"
                        # print("E_ID:", e, " SEQUENCE DO NOT MATCH")
                        # print("API:",iso_seq)
                        # print("MSA:",iso_seq_msa.replace("-",""))
                        total_mismatches += 1
                        if iso_seq == iso_seq_msa.replace("-","")[:-1]:
                            total_mismatches_1 += 1
                    else:
                        total_matches += 1
                    pw2 = pairwise2.align.globalms(wt_seq, iso_seq, 2, -5, -10, -.5)
                    aln_ref = pw2[0][0]
                    aln_isoform = pw2[0][1]
                    if aln_isoform!=iso_seq_msa_corrected:
                        isoforms_with_issue[p+"_"+e] = "Alignment differs than pairwise, see alignment for sanity"
                        total_align_mismatch += 1
                        # print('misalign')
                        # print(aln_isoform)
                        # print(iso_seq_msa_corrected)
                    else:
                        total_align_match += 1

                    gaps = 0
                    gaps_iso = 0
                    missing_pos = []
                    missing_pos_iso = []
                    res_correct = {}
                    isoform_missing_segment = {}
                    count_segment = {}
                    # print("length",len(aln_ref),len(aln_isoform))
                    for i, r in enumerate(aln_ref, 1):
                        if aln_isoform[i-1]=='-':
                            gaps_iso += 1
                        if r == "-":
                            res_correct[i] = ['','','']
                            gaps += 1
                            if aln_isoform[i-1]!="-":
                                # Ref is missing
                                missing_pos_iso.append(i-gaps_iso)
                                if i-gaps==0:
                                    # Take N-term if it's begining
                                    isoform_missing_segment[i-gaps_iso] = r_lookup[1][0]
                                else:
                                    isoform_missing_segment[i-gaps_iso] = (i-gaps,r_lookup[i-gaps][0])
                        else:
                            res_correct[i] = aln_ref[i-gaps-1]
                            if aln_isoform[i-1]=="-":
                                # Ref is missing
                                missing_pos.append(i-gaps)
                            else:
                                segment = r_lookup[i-gaps][0]
                                if segment not in count_segment:
                                    count_segment[segment] = 0
                                count_segment[segment] += 1  

                    result_segment = {}
                    for segment, value in r_segment.items():
                        if segment in count_segment:
                            freq = round(count_segment[segment]/value,2)
                            count = count_segment[segment]
                        else:
                            freq = 0
                            count = 0
                        if freq!=1:
                            # If incomplete segment, save it
                            result_segment[segment] = [freq, count, value]
                    # print(result_segment)
                    # print(missing_pos,missing_pos_iso)
                    ranges = {}
                    ranges['deleted_ref'] = []
                    ranges['inserts'] = []
                    ranges['segments_altered'] = result_segment
                    for k, g in groupby(enumerate(missing_pos), lambda x:x[0]-x[1]):
                        group = list(map(itemgetter(1), g))
                        # What was the previous postions slug
                        from_segment = r_lookup[group[0]][0]
                        to_segment = r_lookup[group[-1]][0]
                        ranges['deleted_ref'].append({'from': [group[0],from_segment] , 'to':[group[-1],to_segment],'length': len(group)})

                    for k, g in groupby(enumerate(missing_pos_iso), lambda x:x[0]-x[1]):
                        group = list(map(itemgetter(1), g))
                        inserted_into = isoform_missing_segment[group[0]]
                        ranges['inserts'].append({'from':group[0], 'to':group[-1], 'inserted_into':inserted_into,'length': len(group)})

                    # print(ranges)

                else:
                    print(e,'no info')

                key = '{}_{}'.format(p,isoform_id)
                dump[key] = ranges
                dump[key]['e_ids'] = es
        #print(dump)
        f = open('protein/data/isoforms.json', 'w')
        json.dump(dump,f, indent=4, separators=(',', ': '))
        print("SUMMARY")
        print("TOTAL MATCHES of isoform seq",total_matches)
        print("ALIGNMENT MATCH",total_align_match,"MISMATCH",total_align_mismatch)
        print("TOTAL MISMATCHES of isoform seq",total_mismatches)
        print("TOTAL MISMATCHES of isoform seq (MSA has one extra)",total_mismatches_1)

        for e,r in isoforms_with_issue.items():
            print(e,r)
    # print(aln_human)
    # print(fasta_raw)
    # data['wt2']=fasta[1]
    # data['pre_aligned']=fasta[1+int(iso)*2]


                            
