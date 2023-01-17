from django.conf import settings
from django.db.models import Q
from django.db import IntegrityError

from protein.models import Protein, ProteinAnomaly, ProteinSegment
from residue.models import Residue, ResidueGenericNumber, ResidueNumberingScheme, ResidueGenericNumberEquivalent

import logging
from collections import OrderedDict
import yaml
import shlex
import os
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

def parse_scheme_tables(path):
    # get generic residue numbering schemes
    rnss = ResidueNumberingScheme.objects.all()
    
    # dictionary for mapping tables
    schemes = OrderedDict()

    for scheme in rnss:
        schemes[scheme.slug] = {}
        schemes[scheme.slug]['generic_numbers'] = {}
        if scheme.parent:
            schemes[scheme.slug]['parent'] = scheme.parent.slug
        schemes[scheme.slug]['obj'] = scheme
        mapping_file = os.sep.join([path, 'mapping_' + scheme.slug + '.txt'])
        if os.path.isfile(mapping_file):
            with open(mapping_file, "r", encoding='UTF-8') as scheme_table_file:
                schemes[scheme.slug]['table'] = {}
                for row in scheme_table_file:
                    split_row = shlex.split(row)
                    schemes[scheme.slug]['table'][split_row[0]] = split_row[1]
    return schemes

def load_reference_positions(path):
    try:
        with open(path, 'r') as ref_position_file:
            ref_positions = yaml.load(ref_position_file, Loader=yaml.FullLoader)
        return ref_positions
    except:
        return False

def create_or_update_residue(protein_conformation, segment, schemes,residue,b_and_c):
    logger = logging.getLogger('build')

    rns_defaults = {'protein_segment': segment} # default numbering scheme for creating generic numbers

    # default numbering scheme
    ns = settings.DEFAULT_NUMBERING_SCHEME
    ns_obj = ResidueNumberingScheme.objects.get(slug=ns)
    
    rvalues = {}
    rvalues['protein_segment'] = segment
    rvalues['amino_acid'] = residue['aa']
    rvalues['generic_number'] = None
    rvalues['display_generic_number'] = None
    sequence_number = residue['pos']
    numbers = residue['numbers']

    if 'generic_number' in numbers:
        numbers = format_generic_numbers(protein_conformation.protein.residue_numbering_scheme, schemes,
                    sequence_number, numbers['generic_number'], numbers['bw'],b_and_c)
        # print(numbers)
    # print(residues,numbers)
    # main generic number
    if 'generic_number' in numbers:
        gnl = numbers['generic_number']
        if gnl in schemes[ns]['generic_numbers']:
            rvalues['generic_number'] = schemes[ns]['generic_numbers'][gnl]
        else:
            try:
                gn, created = ResidueGenericNumber.objects.get_or_create(
                    scheme=ns_obj, label=gnl, defaults=rns_defaults)
                # if created:
                #     logger.info('Created generic number {}'.format(gn.label))
            except IntegrityError:
                gn = ResidueGenericNumber.objects.get(scheme=ns_obj, label=gnl)
            rvalues['generic_number'] = schemes[ns]['generic_numbers'][gnl] = gn


     # equivalent to main generic number
    if 'equivalent' in numbers:
        try:
            gn_equivalent, created = ResidueGenericNumberEquivalent.objects.get_or_create(
                default_generic_number=rvalues['generic_number'],
                scheme=protein_conformation.protein.residue_numbering_scheme,
                defaults={'label': numbers['equivalent']})
            # if created:
            #     logger.info('Created generic number equivalent {} ({}) for scheme {}'.format(
            #         numbers['equivalent'], numbers['generic_number'],
            #         protein_conformation.protein.residue_numbering_scheme))
        except:
            sleep(0.5)
            print('sleep to hope to fix')
            gn_equivalent, created = ResidueGenericNumberEquivalent.objects.get_or_create(
                default_generic_number=rvalues['generic_number'],
                scheme=protein_conformation.protein.residue_numbering_scheme,
                defaults={'label': numbers['equivalent']})
            # gn_equivalent = ResidueGenericNumberEquivalent.objects.get(
            #     default_generic_number=rvalues['generic_number'],
            #     scheme=protein_conformation.protein.residue_numbering_scheme)
    
    # display generic number
    if 'display_generic_number' in numbers:
        ns = protein_conformation.protein.residue_numbering_scheme.slug
        gnl = numbers['display_generic_number']
        if gnl in schemes[ns]['generic_numbers']:
            rvalues['display_generic_number'] = schemes[ns]['generic_numbers'][gnl]
        else:
            try:
                gn, created = ResidueGenericNumber.objects.get_or_create(
                    scheme=protein_conformation.protein.residue_numbering_scheme,
                    label=gnl, defaults=rns_defaults)
                # if created:
                #     logger.info('Created display generic number {}'.format(gn.label))
            except IntegrityError:
                gn = ResidueGenericNumber.objects.get(
                    scheme=protein_conformation.protein.residue_numbering_scheme, label=gnl)
            rvalues['display_generic_number'] = schemes[ns]['generic_numbers'][gnl] = gn
            

        # UPDATE or CREATE the residue
    # bulk_r = Residue(protein_conformation=protein_conformation,sequence_number=sequence_number,defaults = rvalues)
    bulk_r = Residue(protein_conformation=protein_conformation,sequence_number=sequence_number, amino_acid = rvalues['amino_acid'],
        display_generic_number = rvalues['display_generic_number'],
        generic_number = rvalues['generic_number'], protein_segment = segment)
    # r, created = Residue.objects.update_or_create(protein_conformation=protein_conformation,
    #     sequence_number=sequence_number, defaults = rvalues)


    # alternative generic numbers
    # r.alternative_generic_numbers.clear() # remove any existing relations
    bulk_add_alt = []
    if (numbers and 'alternative_generic_numbers' in numbers):
        for alt_scheme, alt_num in numbers['alternative_generic_numbers'].items():
            if alt_num in schemes[alt_scheme]['generic_numbers']:
                argn = schemes[alt_scheme]['generic_numbers'][alt_num]
            else:
                try:
                    argn, created = ResidueGenericNumber.objects.get_or_create(
                        scheme=ResidueNumberingScheme.objects.get(slug=alt_scheme), label=alt_num,
                        defaults=rns_defaults)
                except IntegrityError:
                    argn = ResidueGenericNumber.objects.get(
                        scheme=ResidueNumberingScheme.objects.get(slug=alt_scheme), label=alt_num)
                schemes[alt_scheme]['generic_numbers'][alt_num] = argn
            try:
                bulk_add_alt.append(argn)
                # r.alternative_generic_numbers.add(argn)
            except IntegrityError:
                print('argn already added?')
                pass
                # print('argn already added?')

    return [bulk_r,bulk_add_alt]


def get_gprotein_non_gns(consensus_prot_conf, segment, residues_to_update):
    # Return non-GNs for longest segment in the subfamily
    proteins = Protein.objects.filter(family__parent=consensus_prot_conf.protein.family, species__common_name='Human', accession__isnull=False)
    s_len = 0
    non_gns = []
    for p in proteins:
        resis = Residue.objects.filter(protein_segment=segment, protein_conformation__protein=p)
        if len(resis)>s_len:
            s_len = len(resis)
            non_gns = [r.display_generic_number for r in resis]
    # In case of multiple with same length but different numbers, get all numbers (e.g. Gi/o hbhc)
    if len(residues_to_update)!=len(non_gns):
        new_non_gns = []
        resis = Residue.objects.filter(protein_segment=segment, protein_conformation__protein__species__common_name='Human', 
                                       protein_conformation__protein__accession__isnull=False, protein_conformation__protein__family__parent=consensus_prot_conf.protein.family)
        for r in resis:
            if r.display_generic_number not in new_non_gns:
                new_non_gns.append(r.display_generic_number)
        non_gns = sorted(new_non_gns, key=lambda x:x.label)
    return non_gns


def create_or_update_residues_in_segment(protein_conformation, segment, start, aligned_start, end, aligned_end,
    schemes, ref_positions, protein_anomalies, disregard_db_residues, signprot=False):
    logger = logging.getLogger('build')
    rns_defaults = {'protein_segment': segment} # default numbering scheme for creating generic numbers

    # fetch the residues that should be updated
    residues_to_update = Residue.objects.filter(Q(sequence_number__gte=start) & Q(sequence_number__lte=end),
        protein_conformation=protein_conformation).values_list('sequence_number', 'amino_acid')

    # if the residue records have not been created, use the sequence string instead
    if disregard_db_residues or not residues_to_update.exists():
        residues_to_update = []
        for i, aa in enumerate(protein_conformation.protein.sequence[(start-1):end]):
            residues_to_update.append((start+i, aa)) # replicate the format from values_list

    # how many residues are there before and after the aligned segment (only relevant to non-fully aligned segments)
    if aligned_start:
        residues_before = aligned_start - start
    else:
        residues_before = 0
    if aligned_end:
        residues_after = len(residues_to_update) - (end - aligned_end)
    else:
        residues_after = 0 

    # default numbering scheme
    ns = settings.DEFAULT_NUMBERING_SCHEME
    ns_obj = ResidueNumberingScheme.objects.get(slug=ns)
    
    if signprot:
        non_gns = get_gprotein_non_gns(protein_conformation, segment, residues_to_update)
    created_residues = 0
    for res_num, residue in enumerate(residues_to_update, start=1):
        sequence_number = residue[0]
        
        # dictionary of values for updating/creating a residue
        rvalues = {}
        rvalues['protein_segment'] = segment
        rvalues['amino_acid'] = residue[1]

        # generic numbers
        numbers = None
        rvalues['generic_number'] = None
        rvalues['display_generic_number'] = None

        if (segment.slug in settings.REFERENCE_POSITIONS
            and settings.REFERENCE_POSITIONS[segment.slug] in ref_positions
            and res_num > residues_before 
            and res_num <= residues_after):
            numbers = format_generic_numbers_old(protein_conformation.protein.residue_numbering_scheme, schemes,
                sequence_number, settings.REFERENCE_POSITIONS[segment.slug],
                ref_positions[settings.REFERENCE_POSITIONS[segment.slug]], protein_anomalies)
            
            # main generic number
            if 'generic_number' in numbers:
                gnl = numbers['generic_number']
                if gnl in schemes[ns]['generic_numbers']:
                    rvalues['generic_number'] = schemes[ns]['generic_numbers'][gnl]
                else:
                    try:
                        gn, created = ResidueGenericNumber.objects.get_or_create(
                            scheme=ns_obj, label=gnl, defaults=rns_defaults)
                        if created:
                            logger.info('Created generic number {}'.format(gn.label))
                    except IntegrityError:
                        gn = ResidueGenericNumber.objects.get(scheme=ns_obj, label=gnl)
                    rvalues['generic_number'] = schemes[ns]['generic_numbers'][gnl] = gn

            # equivalent to main generic number
            if 'equivalent' in numbers:
                try:
                    gn_equivalent, created = ResidueGenericNumberEquivalent.objects.get_or_create(
                        default_generic_number=rvalues['generic_number'],
                        scheme=protein_conformation.protein.residue_numbering_scheme,
                        defaults={'label': numbers['equivalent']})
                    if created:
                        logger.info('Created generic number equivalent {} ({}) for scheme {}'.format(
                            numbers['equivalent'], numbers['generic_number'],
                            protein_conformation.protein.residue_numbering_scheme))
                except IntegrityError:
                    try:
                        print("Get",rvalues['generic_number'],protein_conformation.protein.residue_numbering_scheme,)
                        gn_equivalent = ResidueGenericNumberEquivalent.objects.get(
                            default_generic_number=rvalues['generic_number'],
                            scheme=protein_conformation.protein.residue_numbering_scheme)
                    except:
                        print("failed gn_equivalent", rvalues['generic_number'], protein_conformation.protein.residue_numbering_scheme)
            
            # display generic number
            if 'display_generic_number' in numbers:
                ns = protein_conformation.protein.residue_numbering_scheme.slug
                gnl = numbers['display_generic_number']
                if gnl in schemes[ns]['generic_numbers']:
                    rvalues['display_generic_number'] = schemes[ns]['generic_numbers'][gnl]
                else:
                    try:
                        gn, created = ResidueGenericNumber.objects.get_or_create(
                            scheme=protein_conformation.protein.residue_numbering_scheme,
                            label=gnl, defaults=rns_defaults)
                        if created:
                            logger.info('Created generic number {}'.format(gn.label))
                    except IntegrityError:
                        gn = ResidueGenericNumber.objects.get(
                            scheme=protein_conformation.protein.residue_numbering_scheme, label=gnl)
                    rvalues['display_generic_number'] = schemes[ns]['generic_numbers'][gnl] = gn
        elif segment.slug in ProteinSegment.objects.filter(proteinfamily=signprot).values_list('slug', flat=True):
            # if protein_conformation.protein.entry_name!='alpha-consensus':
            rvalues['display_generic_number'] = non_gns[res_num-1]
            rvalues['generic_number'] = non_gns[res_num-1]

        # UPDATE or CREATE the residue
        r, created = Residue.objects.update_or_create(protein_conformation=protein_conformation,
            sequence_number=sequence_number, defaults = rvalues)
        if created:
            created_residues += 1

        # alternative generic numbers
        r.alternative_generic_numbers.clear() # remove any existing relations
        if (segment.slug in settings.REFERENCE_POSITIONS
            and settings.REFERENCE_POSITIONS[segment.slug] in ref_positions
            and numbers and 'alternative_generic_numbers' in numbers):
            for alt_scheme, alt_num in numbers['alternative_generic_numbers'].items():
                if alt_num in schemes[alt_scheme]['generic_numbers']:
                    argn = schemes[alt_scheme]['generic_numbers'][alt_num]
                else:
                    try:
                        argn, created = ResidueGenericNumber.objects.get_or_create(
                            scheme=ResidueNumberingScheme.objects.get(slug=alt_scheme), label=alt_num,
                            defaults=rns_defaults)
                    except IntegrityError:
                        argn = ResidueGenericNumber.objects.get(
                            scheme=ResidueNumberingScheme.objects.get(slug=alt_scheme), label=alt_num)
                    schemes[alt_scheme]['generic_numbers'][alt_num] = argn
                r.alternative_generic_numbers.add(argn)

    if created_residues:
        logger.info('Created {} residues for {} of {}'.format(created_residues, segment, protein_conformation))


def format_generic_numbers_old(residue_numbering_scheme, schemes, sequence_number, ref_position, ref_residue,
    protein_anomalies):
    logger = logging.getLogger('build')

    # empty dict for return values
    numbers = {}

    # generic index
    sgn = ref_position.split("x")
    segment_index = sgn[0]
    ref_generic_index = int(sgn[1])
    generic_index = ref_generic_index - (ref_residue - sequence_number)

    # get generic numbers in display residue numbering scheme
    if ref_position in schemes[residue_numbering_scheme.slug]['table']:
        scheme_ref_position = schemes[residue_numbering_scheme.slug]['table'][ref_position]
        sgn = scheme_ref_position.split("x")
        scheme_ref_generic_index = int(sgn[1])
        scheme_generic_index = scheme_ref_generic_index - (ref_residue - sequence_number)

    # order anomalies if there are more than one
    # this is important for counting offset
    if len(protein_anomalies) > 1:
        protein_anomalies.sort(key=lambda x: x.generic_number.label)
    
    # offset by anomaly (anomalies before and after the reference position are handled differently)
    offset = 0
    prime = ''
    
    # reverse anomalies if the residue position is before the reference position
    # use numbers from the display scheme for this, because the residue may be on the other side of the reference in
    # that scheme
    if scheme_generic_index < ref_generic_index:
        protein_anomalies.reverse()

    for pa in protein_anomalies:
        # generic number without the prime for bulges # FIXME GPCR specific
        spgn = pa.generic_number.label.split("x")
        pa_generic_index = int(spgn[1][:2])
        pa_generic_number = spgn[0] + 'x' + spgn[1][:2]

        # get the anomaly generic number in the display residue numbering scheme
        if pa_generic_number in schemes[residue_numbering_scheme.slug]['table']:
            scheme_pa_generic_number = schemes[residue_numbering_scheme.slug]['table'][pa_generic_number]
            scheme_pa_generic_index = int(scheme_pa_generic_number.split("x")[1])
        
        # add prime to bulges
        if pa.anomaly_type.slug == 'bulge' and (
            (scheme_pa_generic_index < ref_generic_index and scheme_pa_generic_index == (scheme_generic_index+offset)
            or (scheme_pa_generic_index > ref_generic_index
            and scheme_pa_generic_index == (scheme_generic_index+offset-1)))):
            prime = '1'

        # offsets
        if (scheme_pa_generic_index > ref_generic_index and (scheme_generic_index+offset) > scheme_pa_generic_index and
            pa.anomaly_type.slug == 'bulge'):
            offset -= 1
        elif (scheme_pa_generic_index > ref_generic_index and (scheme_generic_index+offset) >= scheme_pa_generic_index
            and pa.anomaly_type.slug == 'constriction'):
            offset += 1
        elif (scheme_pa_generic_index < ref_generic_index and (scheme_generic_index+offset) < scheme_pa_generic_index
            and pa.anomaly_type.slug == 'bulge'):
            offset += 1
        elif (scheme_pa_generic_index < ref_generic_index and (scheme_generic_index+offset) <= scheme_pa_generic_index
            and pa.anomaly_type.slug == 'constriction'):
            offset -= 1

    # structure corrected index (based on anomalies)
    structure_corrected_generic_index = generic_index + offset
    generic_number = segment_index + "x" + str(generic_index)
    structure_corrected_generic_number = segment_index + "x" + str(structure_corrected_generic_index)
    
    # does this generic number exists in the reference tables (relevant to non-fully aligned segments)
    if structure_corrected_generic_number in schemes[residue_numbering_scheme.slug]['table']:
        numbers['generic_number'] = structure_corrected_generic_number + prime
    else:
        return numbers

    # update generic numbers equivalents
    if 'table' in schemes[residue_numbering_scheme.slug]:
        if structure_corrected_generic_number in schemes[residue_numbering_scheme.slug]['table']:
            equivalent = schemes[residue_numbering_scheme.slug]['table'][structure_corrected_generic_number]
            numbers['equivalent'] = equivalent + prime
        else:
            logger.error('{} equivalent for number {} not found, using {}'.format(residue_numbering_scheme.slug,
                structure_corrected_generic_number, structure_corrected_generic_number))

    # alternative schemes
    numbers['alternative_generic_numbers'] = {}
    for scheme, s in schemes.items():
        if 'table' in s:
            if 'parent' in s:
                if generic_number in schemes[s['parent']]['table']:
                    seq_based_num = schemes[s['parent']]['table'][generic_number]
                else:
                    logger.warning('{} equivalent for number {} not found, using {}'.format(s['parent'],
                        generic_number, generic_number))
                    continue
                
                if structure_corrected_generic_number in s['table']:
                    str_based_num = s['table'][structure_corrected_generic_number] + prime
                else:
                    logger.warning('{} equivalent for number {} not found, using {}'.format(scheme,
                        structure_corrected_generic_number, structure_corrected_generic_number))
                    continue
                
                split_str_based_num = str_based_num.split('x')
                numbers['alternative_generic_numbers'][scheme] = seq_based_num + 'x' + split_str_based_num[1]
            elif generic_number in s['table']:
                numbers['alternative_generic_numbers'][scheme] = s['table'][generic_number]

            if scheme == residue_numbering_scheme.slug:
                numbers['display_generic_number'] = numbers['alternative_generic_numbers'][scheme]

    return numbers

def format_generic_numbers(residue_numbering_scheme, schemes, sequence_number, generic_number, bw_number,b_and_c):
    logger = logging.getLogger('build')

    # empty dict for return values
    numbers = {}
    
    # update generic numbers equivalents
    # if 'table' in schemes[residue_numbering_scheme.slug]:
    #     if structure_corrected_generic_number in schemes[residue_numbering_scheme.slug]['table']:
    #         equivalent = schemes[residue_numbering_scheme.slug]['table'][structure_corrected_generic_number]
    #         numbers['equivalent'] = equivalent + prime
    #     else:
    #         logger.error('{} equivalent for number {} not found, using {}'.format(residue_numbering_scheme.slug,
    #             structure_corrected_generic_number, structure_corrected_generic_number))

    # print("GN",generic_number,"bw",bw_number)
    numbers['generic_number'] = generic_number
    generic_number = bw_number.replace(".","x")
    segment_index = generic_number.split("x")[0]
    generic_index = generic_number.split("x")[1]
    if segment_index in b_and_c:
        anomalies = b_and_c[segment_index] #only select segment
    else:
        anomalies = [] #empty

    # alternative schemes
    numbers['alternative_generic_numbers'] = {}
    for scheme, s in schemes.items():
        #print(scheme)
        if 'table' in s:
            if 'parent' in s:
                if generic_number in schemes[s['parent']]['table']:
                    seq_based_num = schemes[s['parent']]['table'][generic_number]
                else:
                    logger.warning('{} equivalent for number {} not found, using {}'.format(s['parent'],
                        generic_number, generic_number))
                    continue

                scheme_anomalies = []
                for anomality in anomalies:
                    pa_generic_number = segment_index + "x" + anomality[:2]
                    #print("pa_generic_number",pa_generic_number)
                    if pa_generic_number in schemes[scheme]['table']:
                        scheme_pa_generic_number = schemes[scheme]['table'][pa_generic_number]
                        scheme_pa_generic_index = int(scheme_pa_generic_number.split("x")[1])
                    if len(anomality)==3: #if it is bulge
                        if int(scheme_pa_generic_index)<50 and int(anomality[:2])>50:
                            #if scheme has bulge before 50 and annotation has it after. Fix offset
                            scheme_pa_generic_index += 1
                        elif int(scheme_pa_generic_index)>50 and int(anomality[:2])<50:
                            #if scheme has bulge after 50 and annotation has it before. Fix offset
                            scheme_pa_generic_index -= 1
                        scheme_pa_generic_index = str(scheme_pa_generic_index) + "1"

                    #print("pa_generic_number",anomality[:2],"scheme_pa_generic_index",scheme_pa_generic_index[:2])
                    #if int(scheme_pa_generic_index[:2])<50 and int(anomality[:2])>50
                    scheme_anomalies.append(str(scheme_pa_generic_index))
                    #print("scheme_pa_generic_index", scheme_pa_generic_index)
                #print(seq_based_num)

                seq_based_num = seq_based_num.replace("x",".") #fix when lookup fails and it uses x based number
                scheme_generic_index = seq_based_num.split(".")[1]
                str_based_num = format_anomalities(scheme_anomalies,scheme_generic_index)
                
                split_str_based_num = str_based_num.split('x')
                numbers['alternative_generic_numbers'][scheme] = seq_based_num + 'x' + str_based_num
            elif generic_number in s['table']:
                numbers['alternative_generic_numbers'][scheme] = s['table'][generic_number]

            if scheme == residue_numbering_scheme.slug:
                numbers['display_generic_number'] = numbers['alternative_generic_numbers'][scheme]
                ### ALSO SET EQUIVALENT
                # equivalent = schemes[residue_numbering_scheme.slug]['table'][structure_corrected_generic_number]
                numbers['equivalent'] = segment_index + "x" + str_based_num

    # if numbers['generic_number'] in ['5x40','5x41','5x42','5x411','5x461','5x43','5x44','5x45','5x46','5x47','5x48','5x49','5x50']:
    #     print(numbers)
    return numbers

def align_protein_to_reference(protein, tpl_ref_pos_file_path, ref_protein):
    logger = logging.getLogger('build')

    # does the template reference position file exists?
    if not os.path.isfile(tpl_ref_pos_file_path):
        logger.error("File {} not found, skipping!".format(tpl_ref_pos_file_path))
        return False
    template_ref_positions = load_reference_positions(tpl_ref_pos_file_path)

    # write sequences to files
    seq_filename = "/tmp/" + protein['entry_name'] + ".fa"
    with open(seq_filename, 'w') as seq_file:
        seq_file.write("> ref\n")
        seq_file.write(ref_protein.sequence + "\n")
        seq_file.write("> seq\n")
        seq_file.write(protein['sequence'] + "\n")

    try:
        ali_filename = "/tmp/out.fa"
        acmd = ClustalOmegaCommandline(infile=seq_filename, outfile=ali_filename, force=True)
        stdout, stderr = acmd()
        a = AlignIO.read(ali_filename, "fasta")
        logger.info("{} aligned to {}".format(protein['entry_name'], ref_protein.entry_name))
    except:
        logger.error('Alignment failed for {}'.format(protein['entry_name']))
        return False

    # find reference positions
    ref_positions = {}
    ref_positions_in_ali = {}
    for position_generic_number, rp in template_ref_positions.items():
        gaps = 0
        for i, r in enumerate(a[0].seq, 1):
            if r == "-":
                gaps += 1
            if i-gaps == rp:
                ref_positions_in_ali[position_generic_number] = i
    for position_generic_number, rp in ref_positions_in_ali.items():
        gaps = 0
        for i, r in enumerate(a[1].seq, 1):
            if r == "-":
                gaps += 1
            if i == rp:
                ref_positions[position_generic_number] = i - gaps

    return ref_positions

def generic_number_within_segment_borders(generic_number, list_of_tpl_generic_numbers):
    generic_index = generic_number.split('x')[1][:2]
    start_index = list_of_tpl_generic_numbers[0].split('x')[1][:2]
    end_index = list_of_tpl_generic_numbers[len(list_of_tpl_generic_numbers)-1].split('x')[1][:2]

    if generic_index >= start_index and generic_index <= end_index:
        return True
    else:
        return False

def format_anomalities(b_and_c,number):
    offset = 0
    bulge = False

    bcs = sorted(b_and_c)
    if int(number)<50:
        bcs = sorted(bcs, reverse=True)
    for bc in bcs:
        if len(bc)>2: #bulge
            # print(bc[0:2],number)
            if int(bc[0:2])<50 and int(number)+offset<int(bc[0:2]): #before x50 and before bulge, do smt
                offset += 1 #eg if 5x461, then 5.46 becomes 5x461, 5.45 becomes 5x46
            elif int(bc[0:2])<50 and int(number)+offset==int(bc[0:2]): #before x50 and is bulge, do smt
                bulge = True # eg if 5x461, then 5.46 becomes 5x461
            elif int(bc[0:2])>=50 and int(number)+offset>int(bc[0:2])+1: #after x50 and after bulge, do smt
                offset -= 1 #eg if 2x551, then 2.56 becomes 2x551, 5.57 becomes 5x56
            elif int(bc[0:2])>=50 and int(number)+offset==int(bc[0:2])+1: #after x50 and 1 after bulge, do smt
                bulge = True # eg if 2x551, then 2.56 becomes 2x551

        else: #2 numbers, it's a constriction
            if int(bc[0:2])<50 and int(number)+offset<=int(bc[0:2]): #before x50 and before or equal constrictions, do smt
                offset -= 1 #eg if constriction is 7x44, then 7.44 becomes 7x43, 7.43 becomes 7x42
            if int(bc[0:2])>50 and int(number)+offset>=int(bc[0:2]): #before x50 and before or equal constrictions, do smt
                offset += 1 #eg if constriction is 4x57, then 4.57 becomes 4x58, 4.58 becomes 4x59

    if bulge!=True:
        gn = str(int(number)+offset)
    elif int(number)<50:
        gn = str(int(number)+offset)+"1"
    elif int(number)>=50:
        gn = str(int(number)-1+offset)+"1"
    elif int(number)==50:
        print('BULGE IN x50',number,offset)
        gn = str(int(number))+"1"

    return gn

def ggn(gn):
    ''' Converts display generic number to generic number.
    '''
    return gn.split('.')[0]+'x'+gn.split('x')[1]

def dgn(gn, protein_conformation):
    ''' Converts generic number to display generic number.
    '''
    scheme = ResidueNumberingScheme.objects.get(slug=protein_conformation.protein.residue_numbering_scheme.slug)
    convert_gn = ResidueGenericNumberEquivalent.objects.get(label=gn, scheme=scheme).default_generic_number.label
    return Residue.objects.get(protein_conformation=protein_conformation, generic_number__label=convert_gn).display_generic_number.label
    
class DummyResidue():
    def __init__(self, sequence_number, amino_acid, three_letter, display_generic_number=None, generic_number=None, protein_conformation=None, protein_segment=None):
        self.sequence_number = sequence_number
        self.amino_acid = amino_acid
        self.three_letter = three_letter
        self.display_generic_number = display_generic_number
        self.generic_number = generic_number
        self.protein_conformation = protein_conformation
        self.protein_segment = protein_segment

    def __repr__(self):
        return self.amino_acid+str(self.sequence_number)
