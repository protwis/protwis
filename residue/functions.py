from django.conf import settings

from protein.models import ProteinAnomaly
from residue.models import Residue, ResidueGenericNumber, ResidueNumberingScheme

import logging
from collections import OrderedDict
import yaml
import shlex
import os

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
    with open(path, 'r') as ref_position_file:
        ref_positions = yaml.load(ref_position_file)
    return ref_positions

def create_or_update_residues_in_segment(protein_conformation, segment, start, end, schemes, ref_positions,
    protein_anomalies):
    logger = logging.getLogger(__name__)
    rns_defaults = {'protein_segment': segment} # default numbering scheme for creating generic numbers
    for i, aa in enumerate(protein_conformation.protein.sequence[(start-1):end]):
        sequence_number = start + i
        
        # dictionary of values for updating/creating a residue
        rvalues = {}
        rvalues['protein_segment'] = segment
        rvalues['amino_acid'] = aa    

        # generic numbers
        if (segment.slug in settings.REFERENCE_POSITIONS and
            settings.REFERENCE_POSITIONS[segment.slug] in ref_positions):
            numbers = format_generic_numbers(protein_conformation.protein.residue_numbering_scheme, schemes,
                sequence_number, settings.REFERENCE_POSITIONS[segment.slug],
                ref_positions[settings.REFERENCE_POSITIONS[segment.slug]], protein_anomalies)
            
            # main generic number
            ns = settings.DEFAULT_NUMBERING_SCHEME
            gnl = numbers['generic_number']
            if gnl in schemes[ns]['generic_numbers']:
                rvalues['generic_number'] = schemes[ns]['generic_numbers'][gnl]
            else:
                gn, created = ResidueGenericNumber.objects.get_or_create(
                    scheme=ResidueNumberingScheme.objects.get(slug=ns),
                    label=gnl, defaults=rns_defaults)
                rvalues['generic_number'] = schemes[ns]['generic_numbers'][gnl] = gn
                if created:
                    logger.info('Created generic number {}'.format(gn.label))
            
            # display generic number
            ns = protein_conformation.protein.residue_numbering_scheme.slug
            gnl = numbers['display_generic_number']
            if gnl in schemes[ns]['generic_numbers']:
                rvalues['display_generic_number'] = schemes[ns]['generic_numbers'][gnl]
            else:
                gn, created = ResidueGenericNumber.objects.get_or_create(
                    scheme=protein_conformation.protein.residue_numbering_scheme,
                    label=gnl, defaults=rns_defaults)
                rvalues['display_generic_number'] = schemes[ns]['generic_numbers'][gnl] = gn
                if created:
                    logger.info('Created generic number {}'.format(gn.label))
            
        # UPDATE or CREATE the residue
        r, created = Residue.objects.update_or_create(
            protein_conformation=protein_conformation,
            sequence_number=sequence_number,
            defaults = rvalues)
        if created:
            if r.generic_number:
                logger.info('Created residue {}{}({}) for protein {}'.format(r.amino_acid, r.sequence_number,
                    r.generic_number.label, protein_conformation.protein.entry_name))
            else:
                logger.info('Created residue {}{} for protein {}'.format(r.amino_acid, r.sequence_number,
                    protein_conformation.protein.entry_name))

        # alternative generic numbers
        if (segment.slug in settings.REFERENCE_POSITIONS and
            settings.REFERENCE_POSITIONS[segment.slug] in ref_positions):
            r.alternative_generic_numbers.clear() # remove any existing relations
            for alt_scheme, alt_num in numbers['alternative_generic_numbers'].items():
                if alt_num in schemes[alt_scheme]['generic_numbers']:
                    argn = schemes[alt_scheme]['generic_numbers'][alt_num]
                else:
                    argn, created = ResidueGenericNumber.objects.get_or_create(
                        scheme=ResidueNumberingScheme.objects.get(slug=alt_scheme), label=alt_num,
                        defaults=rns_defaults)
                    schemes[alt_scheme]['generic_numbers'][alt_num] = argn
                r.alternative_generic_numbers.add(argn)

def format_generic_numbers(residue_numbering_scheme, schemes, sequence_number, ref_position, ref_residue,
    protein_anomalies):
    numbers = {}

    # generic index
    sgn = ref_position.split("x")
    ref_generic_index = int(sgn[1])
    generic_index = ref_generic_index - (ref_residue - sequence_number)

    # offset by anomaly (anomalies before and after the reference position are handled differently)
    if len(protein_anomalies) > 1:
        protein_anomalies.sort()
    offset = 0
    prime = ''
    
    # reverse anomalies if the residue position is before the reference position
    if generic_index < ref_generic_index:
        protein_anomalies.reverse()

    for pal in protein_anomalies:
        pa = ProteinAnomaly.objects.get(generic_number__label=pal) # FIXME prefetch
        spagn = pa.generic_number.label.split("x")
        pa_generic_index = int(spagn[1][:2]) # generic number without the prime for bulges
        
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
    
    generic_number = sgn[0] + "x" + str(generic_index)
    structure_corrected_generic_number = sgn[0] + "x" + str(structure_corrected_generic_index)
    numbers['generic_number'] = structure_corrected_generic_number + prime

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