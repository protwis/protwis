from django.conf import settings
from django.db.models import Q

from protein.models import ProteinAnomaly
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
            ref_positions = yaml.load(ref_position_file)
        return ref_positions
    except:
        # if the file does not exists, create an empty file
        ref_position_file = open(path, 'w')
        ref_position_file.close()

def create_or_update_residues_in_segment(protein_conformation, segment, start, end, schemes, ref_positions,
    protein_anomalies, disregard_db_residues):
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

    # default numbering scheme
    ns = settings.DEFAULT_NUMBERING_SCHEME
    ns_obj = ResidueNumberingScheme.objects.get(slug=ns)
    
    updated_residues = 0
    for residue in residues_to_update:
        sequence_number = residue[0]
        
        # dictionary of values for updating/creating a residue
        rvalues = {}
        rvalues['protein_segment'] = segment
        rvalues['amino_acid'] = residue[1]

        # generic numbers
        if (segment.slug in settings.REFERENCE_POSITIONS and
            settings.REFERENCE_POSITIONS[segment.slug] in ref_positions):
            numbers = format_generic_numbers(protein_conformation.protein.residue_numbering_scheme, schemes,
                sequence_number, settings.REFERENCE_POSITIONS[segment.slug],
                ref_positions[settings.REFERENCE_POSITIONS[segment.slug]], protein_anomalies)
            
            # main generic number
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
                    gn_equivalent = ResidueGenericNumberEquivalent.objects.get(
                        default_generic_number=rvalues['generic_number'],
                        scheme=protein_conformation.protein.residue_numbering_scheme)
            
            # display generic number
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
        else:
            rvalues['generic_number'] = None
            rvalues['display_generic_number'] = None
            
        # UPDATE or CREATE the residue
        r, created = Residue.objects.update_or_create(
            protein_conformation=protein_conformation,
            sequence_number=sequence_number,
            defaults = rvalues)
        if created:
            updated_residues += 1

        # alternative generic numbers
        r.alternative_generic_numbers.clear() # remove any existing relations
        if (segment.slug in settings.REFERENCE_POSITIONS and
            settings.REFERENCE_POSITIONS[segment.slug] in ref_positions):
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

    logger.info('Created {} residues for {} of {}'.format(updated_residues, segment, protein_conformation))

def format_generic_numbers(residue_numbering_scheme, schemes, sequence_number, ref_position, ref_residue,
    protein_anomalies):
    numbers = {}

    # generic index
    sgn = ref_position.split("x")
    segment_index = sgn[0]
    ref_generic_index = int(sgn[1])
    generic_index = ref_generic_index - (ref_residue - sequence_number)

    # order anomalies if there are more than one
    # this is important for counting offset
    if len(protein_anomalies) > 1:
        protein_anomalies.sort(key=lambda x: x.generic_number.label)
    
    # offset by anomaly (anomalies before and after the reference position are handled differently)
    offset = 0
    prime = ''
    
    # reverse anomalies if the residue position is before the reference position
    if generic_index < ref_generic_index:
        protein_anomalies.reverse()

    for pa in protein_anomalies:
        pa_generic_index = int(pa.generic_number.label.split("x")[1][:2]) # generic number without the prime for bulges
        
        # add prime to bulges
        if pa.anomaly_type.slug == 'bulge' and (
            (pa_generic_index < ref_generic_index and pa_generic_index == (generic_index+offset) or
            (pa_generic_index > ref_generic_index and pa_generic_index == (generic_index+offset-1)))):
            prime = '1'

        # offsets
        if (pa_generic_index > ref_generic_index and (generic_index+offset) > pa_generic_index and
            pa.anomaly_type.slug == 'bulge'):
            offset -= 1
        elif (pa_generic_index > ref_generic_index and (generic_index+offset) >= pa_generic_index and
            pa.anomaly_type.slug == 'constriction'):
            offset += 1
        elif (pa_generic_index < ref_generic_index and (generic_index+offset) < pa_generic_index and
            pa.anomaly_type.slug == 'bulge'):
            offset += 1
        elif (pa_generic_index < ref_generic_index and (generic_index+offset) <= pa_generic_index and
            pa.anomaly_type.slug == 'constriction'):
            offset -= 1

    # structure corrected index (based on anomalies)
    structure_corrected_generic_index = generic_index + offset
    
    generic_number = segment_index + "x" + str(generic_index)
    structure_corrected_generic_number = segment_index + "x" + str(structure_corrected_generic_index)
    numbers['generic_number'] = structure_corrected_generic_number + prime

    # update generic numbers equivalents
    if 'table' in schemes[residue_numbering_scheme.slug]:
        equivalent = schemes[residue_numbering_scheme.slug]['table'][structure_corrected_generic_number]
        numbers['equivalent'] = equivalent + prime

    # alternative schemes
    numbers['alternative_generic_numbers'] = {}
    for scheme, s in schemes.items():
        if 'table' in s:
            if 'parent' in s:
                seq_based_num = schemes[s['parent']]['table'][generic_number]
                str_based_num = s['table'][structure_corrected_generic_number] + prime
                split_str_based_num = str_based_num.split('x')
                numbers['alternative_generic_numbers'][scheme] = seq_based_num + 'x' + split_str_based_num[1]
            elif generic_number in s['table']:
                numbers['alternative_generic_numbers'][scheme] = s['table'][generic_number]

            if scheme == residue_numbering_scheme.slug:
                numbers['display_generic_number'] = numbers['alternative_generic_numbers'][scheme]

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