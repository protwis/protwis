from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django import forms
from django.db.models import Count, Min, Sum, Avg, Q
from django.utils.text import slugify
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from interaction.models import *
from interaction.forms import PDBform
from ligand.models import Ligand
from ligand.models import LigandType
from ligand.models import LigandRole, LigandProperities
from structure.models import Structure, PdbData, Rotamer, Fragment
from structure.functions import BlastSearch
from structure.assign_generic_numbers_gpcr import GenericNumbering
from protein.models import ProteinConformation, Protein, ProteinSegment
from residue.models import Residue, ResidueGenericNumber, ResidueGenericNumberEquivalent, ResidueNumberingScheme
from common.models import WebResource
from common.models import WebLink
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot
from common.selection import SimpleSelection, Selection, SelectionItem
from common import definitions
from common.views import AbsTargetSelection
from common.alignment import Alignment
from protein.models import Protein, ProteinFamily, ProteinGProtein, ProteinGProteinPair

import os
from os import listdir, devnull, makedirs
import yaml
from operator import itemgetter
from datetime import datetime
import re
import json
import logging
from subprocess import call, Popen, DEVNULL
import urllib
import collections
from collections import OrderedDict
from io import StringIO, BytesIO
from Bio.PDB import PDBIO, PDBParser
import xlsxwriter

######@
import numpy as np

AA = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
      'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
      'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
      'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
      'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def regexaa(aa):
    aaPattern = re.compile(r'^(\w{3})(\d+)([\w\s]+)$')
    # Splits the string into AA number CHAIN : LEU339A => ('LEU', '339', 'A')
    result = aaPattern.search(aa)
    if result:
        result = result.groups()
        aa = AA[result[0]]
        number = result[1]
        chain = result[2]
        return aa, number, chain
    else:
        aaPattern = re.compile(r'^(\w{3})(\d+)$')
        result = aaPattern.search(aa)
        if result:
            result = result.groups()
            aa = AA[result[0]]
            number = result[1]
            chain = ' '
            return aa, number, chain
        else:
            return None, None, None


class InteractionSelection(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 1
    docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Select the structure of interest by using the dropdown in the middle. The selection if viewed to the right and the interactions will be loaded immediately.'

    # Middle section
    numbering_schemes = False
    filters = False
    search = False
    title = "Select a structure based on PDB-code"

    template_name = 'interaction/interactionselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # Buttons
    buttons = {
        'continue': {
            'label': 'Show interactions',
            'onclick': 'submitupload()',
            'color': 'success',
        }
    }

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        context['structures'] = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index', 'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate(
            num_ligands=Count('structure_ligand_pair', distinct=True), num_interactions=Count('pk', distinct=True)).order_by('structure_ligand_pair__structure__pdb_code__index')
        context['structure_groups'] = sorted(set([ structure['structure_ligand_pair__structure__pdb_code__index'][0] for structure in context['structures'] ]))
        context['form'] = PDBform()
        return context


def StructureDetails(request, pdbname):
    """
    Show structure details
    """
    pdbname = pdbname
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name', 'structure_ligand_pair__pdb_reference', 'structure_ligand_pair__annotated').filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
    resn_list = ''


    main_ligand = "none"
    for structure in structures:
        if structure['structure_ligand_pair__annotated']:
            resn_list += ",\"" + \
                structure['structure_ligand_pair__pdb_reference'] + "\""
            main_ligand = structure['structure_ligand_pair__pdb_reference']



    crystal = Structure.objects.get(pdb_code__index=pdbname)
    p = Protein.objects.get(protein=crystal.protein_conformation.protein)

    residuelist = Residue.objects.filter(protein_conformation__protein=p).prefetch_related('protein_segment','display_generic_number','generic_number')
    lookup = {}

    residues_lookup = {}
    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number.label] = r.sequence_number
            residues_lookup[r.sequence_number] = r.amino_acid +str(r.sequence_number)+ " "+ r.generic_number.label

    residues = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname, structure_ligand_pair__annotated=True).exclude(interaction_type__type ='hidden').order_by('rotamer__residue__sequence_number')
    residues_browser = []
    ligands = []
    display_res = []
    main_ligand_full = "None"
    residue_table_list = []
    for residue in residues:
        key = residue.interaction_type.name
        aa = residue.rotamer.residue.amino_acid
        pos = residue.rotamer.residue.sequence_number
        wt_pos = -1

        if residue.rotamer.residue.generic_number and residue.rotamer.residue.generic_number.label in lookup:
            residue_table_list.append(
                residue.rotamer.residue.generic_number.label)
            wt_pos = lookup[residue.rotamer.residue.generic_number.label]

        if residue.rotamer.residue.protein_segment:
            segment = residue.rotamer.residue.protein_segment.slug
        else:
            segment = ''
        if residue.rotamer.residue.display_generic_number:
            display = residue.rotamer.residue.display_generic_number.label
        else:
            display = ''
        ligand = residue.structure_ligand_pair.ligand.name
        display_res.append(str(pos))
        residues_browser.append({'type': key, 'aa': aa, 'ligand': ligand,
                                 'pos': pos, 'wt_pos': wt_pos, 'gpcrdb': display, 'segment': segment})

        if pos not in residues_lookup:
            residues_lookup[pos] = aa + str(pos) + " " +display + " interaction " + key
        else:
            residues_lookup[pos] += " interaction " + key

        if ligand not in ligands:
            ligands.append(ligand)
            main_ligand_full = ligand
    display_res = ' or '.join(display_res)
    # RESIDUE TABLE
    segments = ProteinSegment.objects.all().filter().prefetch_related()

    proteins = [p]

    numbering_schemes_selection = [settings.DEFAULT_NUMBERING_SCHEME]
    numbering_schemes_selection.append(p.residue_numbering_scheme.slug)

    numbering_schemes = ResidueNumberingScheme.objects.filter(
        slug__in=numbering_schemes_selection).all()
    default_scheme = numbering_schemes.get(
        slug=settings.DEFAULT_NUMBERING_SCHEME)
    data = OrderedDict()

    for segment in segments:
        data[segment.slug] = OrderedDict()
        residues = Residue.objects.filter(protein_segment=segment,  protein_conformation__protein=p,
                                          generic_number__label__in=residue_table_list).prefetch_related('protein_conformation__protein',
                                                                                                         'protein_conformation__state', 'protein_segment',
                                                                                                         'generic_number', 'display_generic_number', 'generic_number__scheme',
                                                                                                         'alternative_generic_numbers__scheme')
        for scheme in numbering_schemes:
            if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {
                        scheme.slug: pos, 'seq': ['-'] * len(proteins)}
            elif scheme == default_scheme:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {
                        scheme.slug: pos, 'seq': ['-'] * len(proteins)}

        for residue in residues:
            alternatives = residue.alternative_generic_numbers.all()
            pos = residue.generic_number
            for alternative in alternatives:
                if alternative.scheme not in numbering_schemes:
                    continue
                scheme = alternative.scheme
                if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                    pos = residue.generic_number
                    if scheme == pos.scheme:
                        data[segment.slug][pos.label]['seq'][proteins.index(
                            residue.protein_conformation.protein)] = str(residue)
                    else:
                        if scheme.slug not in data[segment.slug][pos.label].keys():
                            data[segment.slug][pos.label][
                                scheme.slug] = alternative.label
                        if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                            data[segment.slug][pos.label][
                                scheme.slug] += " " + alternative.label
                        data[segment.slug][pos.label]['seq'][proteins.index(
                            residue.protein_conformation.protein)] = str(residue)
                else:
                    if scheme.slug not in data[segment.slug][pos.label].keys():
                        data[segment.slug][pos.label][
                            scheme.slug] = alternative.label
                    if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                        data[segment.slug][pos.label][
                            scheme.slug] += " " + alternative.label
                    data[segment.slug][pos.label]['seq'][proteins.index(
                        residue.protein_conformation.protein)] = str(residue)

    # Preparing the dictionary of list of lists. Dealing with tripple nested
    # dictionary in django templates is a nightmare
    flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])
    for s in iter(flattened_data):
        flattened_data[s] = [[data[s][x][y.slug]
                              for y in numbering_schemes] + data[s][x]['seq'] for x in sorted(data[s])]

    context = {}
    context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [
                            x.name for x in proteins], [x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
    context['segments'] = [x.slug for x in segments if len(data[x.slug])]
    context['data'] = flattened_data
    context['number_of_schemes'] = len(numbering_schemes)

    residuelist = Residue.objects.filter(protein_conformation__protein=p).prefetch_related('protein_segment','display_generic_number','generic_number')
    HelixBox = DrawHelixBox(
                residuelist, p.get_protein_class(), str(p), nobuttons=1)
    SnakePlot = DrawSnakePlot(
                residuelist, p.get_protein_class(), str(p), nobuttons=1)

    return render(request, 'interaction/structure.html', {'pdbname': pdbname, 'structures': structures,
                                                          'crystal': crystal, 'protein': p, 'helixbox' : HelixBox, 'snakeplot': SnakePlot, 'residues': residues_browser, 'residues_lookup': residues_lookup, 'display_res': display_res, 'annotated_resn':
                                                          resn_list, 'ligands': ligands,'main_ligand' : main_ligand,'main_ligand_full' : main_ligand_full, 'data': context['data'],
                                                          'header': context['header'], 'segments': context['segments'],
                                                          'number_of_schemes': len(numbering_schemes)})


def list_structures(request):
    form = PDBform()
    #structures = ResidueFragmentInteraction.objects.distinct('structure_ligand_pair__structure').all()
    structures = ResidueFragmentInteraction.objects.exclude(interaction_type__slug='acc').values('structure_ligand_pair__structure__pdb_code__index', 'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate(
        num_ligands=Count('structure_ligand_pair', distinct=True), num_interactions=Count('pk', distinct=True)).order_by('structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name')
    #structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index','structure_ligand_pair__structure').annotate(Count('structure_ligand_pair__ligand'))
    # print(structures.count())
    genes = {}
    countligands = {}
    totalligands = 0
    totalinteractions = 0
    totaltopinteractions = 0
    for structure in structures:
        # print(structure)
        if structure['structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name'] not in genes:
            genes[structure[
                'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name']] = 1
        totalligands += structure['num_ligands']
        totalinteractions += structure['num_interactions']
        ligands = ResidueFragmentInteraction.objects.exclude(interaction_type__slug='acc').values('structure_ligand_pair__ligand__name').filter(structure_ligand_pair__structure__pdb_code__index=structure[
            'structure_ligand_pair__structure__pdb_code__index']).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
        for ligand in ligands:
            totaltopinteractions += ligand['numRes']
            if ligand['structure_ligand_pair__ligand__name'] not in countligands:
                countligands[ligand['structure_ligand_pair__ligand__name']] = 1
            break

        # print(structure.structure_ligand_pair.structure.pdb_code.index)
        # print(structure.numRes)
    #objects = Model.objects.filter(id__in=object_ids)
    #context = {}
    print('Structures with ligand information:' + str(structures.count()))
    print('Distinct genes:' + str(len(genes)))
    #print('ligands:' + str(totalligands))
    print('interactions:' + str(totalinteractions))
    print('interactions from top ligands:' + str(totaltopinteractions))
    print('Distinct ligands as top ligand:' + str(len(countligands)))

    return render(request, 'interaction/list.html', {'form': form, 'structures': structures})


def crystal(request):
    pdbname = request.GET.get('pdb')
    form = PDBform()
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name').filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
    crystal = Structure.objects.get(pdb_code__index=pdbname)
    p = Protein.objects.get(protein=crystal.protein_conformation.protein)
    residues = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).order_by('rotamer__residue__sequence_number')
    print("residues", residues)
    return render(request, 'interaction/crystal.html', {'form': form, 'pdbname': pdbname, 'structures': structures, 'crystal': crystal, 'protein': p, 'residues': residues})


def view(request):
    pdbname = request.GET.get('pdb')
    form = PDBform()
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name').filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
    return render(request, 'interaction/view.html', {'form': form, 'pdbname': pdbname, 'structures': structures})


def ligand(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    form = PDBform()
    fragments = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname).filter(
        structure_ligand_pair__ligand__name=ligand).order_by('interaction_type')
    return render(request, 'interaction/ligand.html', {'form': form, 'pdbname': pdbname, 'ligand': ligand, 'fragments': fragments})


def fragment(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    fragment = request.GET.get('fragment')
    form = PDBform()
    fragments = ResidueFragmentInteraction.objects.get(id=fragment)
    return render(request, 'interaction/fragment.html', {'form': form, 'pdbname': pdbname, 'ligand': ligand, 'fragmentid': fragment, 'fragments': fragments})


def updateall(request):
    structures = Structure.objects.values('pdb_code__index').distinct()
    for s in structures:
        pdbname = s['pdb_code__index']
        check = ResidueFragmentInteraction.objects.filter(
            structure_ligand_pair__structure__pdb_code__index=pdbname).all()

        if check.count() == 0:
            t1 = datetime.now()
            runcalculation(pdbname)
            t2 = datetime.now()
            delta = t2 - t1
            seconds = delta.total_seconds()
            print("Calculation: Total time " +
                  str(seconds) + " seconds for " + pdbname)
            t1 = datetime.now()
            results = parsecalculation(pdbname, False)
            t2 = datetime.now()
            delta = t2 - t1
            seconds = delta.total_seconds()
            print("Parsing: Total time " +
                  str(seconds) + " seconds for " + pdbname)
            check = ResidueFragmentInteraction.objects.filter(
                structure_ligand_pair__structure__pdb_code__index=pdbname).all()
            print("Interactions found: " + str(check.count()))
        else:
            print(pdbname + " already calculated")

    # return render(request,'interaction/view.html',{'form': form, 'pdbname':
    # pdbname, 'structures': structures})


def runcalculation(pdbname, peptide=""):
    calc_script = os.sep.join(
        [os.path.dirname(__file__), 'legacy_functions.py'])

    call(["python2.7", calc_script, "-p", pdbname, "-c", peptide],
         stdout=open(devnull, 'wb'), stderr=open(devnull, 'wb'))

    return None


def check_residue(protein, pos, aa):
    residue = Residue.objects.filter(
        protein_conformation=protein, sequence_number=pos)
    if residue.exists():
        residue = Residue.objects.get(
            protein_conformation=protein, sequence_number=pos)
        if residue.amino_acid != aa:
            residue.amino_acid = aa
            residue.save()
    else:
        residue, created = Residue.objects.get_or_create(
            protein_conformation=protein, sequence_number=pos, amino_acid=aa)
        # continue #SKIP THESE -- mostly fusion residues that aren't mapped
        # yet.
    return residue


def extract_fragment_rotamer(f, residue, structure, ligand):
    if os.path.isfile(f):
        f_in = open(f, 'r')
        rotamer_pdb = ''
        fragment_pdb = ''
        for line in f_in:
            if line.startswith('HETATM') or line.startswith('CONECT') or line.startswith('MASTER') or line.startswith('END'):
                fragment_pdb += line
            elif line.startswith('ATOM'):
                rotamer_pdb += line
            else:
                fragment_pdb += line
                rotamer_pdb += line
        f_in.close()

        rotamer_data, created = PdbData.objects.get_or_create(pdb=rotamer_pdb)
        rotamer, created = Rotamer.objects.get_or_create(
            residue=residue, structure=structure, pdbdata=rotamer_data)

        fragment_data, created = PdbData.objects.get_or_create(
            pdb=fragment_pdb)
        fragment, created = Fragment.objects.get_or_create(
            ligand=ligand, structure=structure, pdbdata=fragment_data, residue=residue)
    else:
        #quit("Could not find " + residue)
        return None, None

    return fragment, rotamer


# consider skipping non hetsym ligands FIXME
def parsecalculation(pdbname, debug=True, ignore_ligand_preset=False):
    logger = logging.getLogger('build')
    mypath = '/tmp/interactions/results/' + pdbname + '/output'
    module_dir = '/tmp/interactions'
    results = []
    web_resource = web_resource = WebResource.objects.get(slug='pdb')
    web_link, created = WebLink.objects.get_or_create(
        web_resource=web_resource, index=pdbname)

    annotated_found = 0

    structure = Structure.objects.filter(pdb_code=web_link)
    if structure.exists():
        structure = Structure.objects.get(pdb_code=web_link)

        if structure.pdb_data is None:
            f = module_dir + "/pdbs/" + pdbname + ".pdb"
            if os.path.isfile(f):
                pdbdata, created = PdbData.objects.get_or_create(
                    pdb=open(f, 'r').read())  # does this close the file?
            else:
                print('quitting due to no pdb in filesystem')
                quit()
            structure.pdb_data = pdbdata
            structure.save()

        protein = structure.protein_conformation

        for f in listdir(mypath):
            if os.path.isfile(os.path.join(mypath, f)):
                annotated = 0
                #print(mypath + "/" +f)
                result = yaml.load(open(mypath + "/" + f, 'rb'), Loader=yaml.FullLoader)
                output = result
                temp = f.replace('.yaml', '').split("_")
                temp.append([output])
                temp.append(round(output['score']))
                temp.append((output['inchikey']).strip())
                temp.append((output['smiles']).strip())

                results.append(temp)

                if 'prettyname' not in output:
                    # use hetsyn name if possible, others 3letter
                    output['prettyname'] = temp[1]

                f = module_dir + "/results/" + pdbname + "/interaction" + \
                    "/" + pdbname + "_" + temp[1] + ".pdb"
                if os.path.isfile(f):
                    pdbdata, created = PdbData.objects.get_or_create(
                        pdb=open(f, 'r').read())  # does this close the file?
                    if debug:
                        print("Found file" + f)
                else:
                    print('quitting due to no pdb for fragment in filesystem', f)
                    quit()

                structureligandinteraction = StructureLigandInteraction.objects.filter(
                    pdb_reference=temp[1], structure=structure, annotated=True) #, pdb_file=None
                if structureligandinteraction.exists():  # if the annotated exists
                    annotated_found = 1
                    annotated = 1
                    try:
                        structureligandinteraction = structureligandinteraction.get()
                        structureligandinteraction.pdb_file = pdbdata
                        ligand = structureligandinteraction.ligand
                        if structureligandinteraction.ligand.properities.inchikey is None:
                            structureligandinteraction.ligand.properities.inchikey = output['inchikey'].strip()
                        elif structureligandinteraction.ligand.properities.inchikey != output['inchikey'].strip():
                            logger.error(
                                'Ligand/PDB inchikey mismatch (PDB:' + pdbname + ' LIG:' + output['prettyname'] + '): '+structureligandinteraction.ligand.properities.inchikey+' vs '+ output['inchikey'].strip())
                    except Exception as msg:
                        print('error with dublication structureligand',temp[1],msg)
                        break
                elif StructureLigandInteraction.objects.filter(pdb_reference=temp[1], structure=structure).exists():
                    try:
                        structureligandinteraction = StructureLigandInteraction.objects.filter(
                            pdb_reference=temp[1], structure=structure).get()
                        structureligandinteraction.pdb_file = pdbdata
                    except: #already there
                        structureligandinteraction = StructureLigandInteraction.objects.filter(
                            pdb_reference=temp[1], structure=structure, pdb_file=pdbdata).get()
                    ligand = structureligandinteraction.ligand
                else:  # create ligand and pair

                    ligand = Ligand.objects.filter(
                        name=output['prettyname'], canonical=True)

                    if ligand.exists():  # if ligand with name (either hetsyn or 3 letter) exists use that.
                        ligand = ligand.get()
                    else:  # create it
                        default_ligand_type = 'N/A'
                        lt, created = LigandType.objects.get_or_create(slug=slugify(default_ligand_type),
                                                                       defaults={'name': default_ligand_type})

                        ligand = Ligand()
                        ligand = ligand.load_from_pubchem(
                            'inchikey', output['inchikey'].strip(), lt, output['prettyname'])
                        try:
                            ligand.save()
                        except:
                            #print('ligand save failed, empty ligand?',output['prettyname'])
                            continue

                    ligandrole, created = LigandRole.objects.get_or_create(
                        name='unknown', slug='unknown')
                    structureligandinteraction = StructureLigandInteraction()
                    structureligandinteraction.ligand = ligand
                    structureligandinteraction.structure = structure
                    structureligandinteraction.ligand_role = ligandrole
                    structureligandinteraction.pdb_file = pdbdata
                    structureligandinteraction.pdb_reference = temp[1]

                structureligandinteraction.save()

                ResidueFragmentInteraction.objects.filter(structure_ligand_pair=structureligandinteraction).delete()

                for interaction in output['interactions']:
                    # print(interaction)
                    aa = interaction[0]
                    aa, pos, chain = regexaa(aa)
                    residue = check_residue(protein, pos, aa)
                    f = interaction[1]

                    fragment, rotamer = extract_fragment_rotamer(
                                    f, residue, structure, ligand)

                    # print(interaction[2],interaction[3],interaction[4],interaction[5])
                    if fragment!=None:
                        interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(
                                        slug=interaction[2], name=interaction[3], type=interaction[4], direction=interaction[5])
                        fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(
                                        structure_ligand_pair=structureligandinteraction, interaction_type=interaction_type, fragment=fragment, rotamer=rotamer)
                #print("Inserted",len(output['interactions']),"interactions","ligand",temp[1],"annotated",annotated)
        # if not annotated_found:
        #     print("No interactions for annotated ligand")

    else:
        if debug:
            logger.info("Structure not in DB?!??!")
        for f in listdir(mypath):
            if os.path.isfile(os.path.join(mypath, f)):
                result = yaml.load(open(mypath + "/" + f, 'rb'), Loader=yaml.FullLoader)
                output = result

                temp = f.replace('.yaml', '').split("_")
                temp.append([output])
                temp.append(round(output['score']))
                temp.append((output['inchikey']).strip())
                temp.append((output['smiles']).strip())
                results.append(temp)

    results = sorted(results, key=itemgetter(3), reverse=True)

    return results


def runusercalculation(filename, session):
    calc_script = os.sep.join(
        [os.path.dirname(__file__), 'legacy_functions.py'])
    call(["python2.7", calc_script, "-p", filename, "-s", session])
    return None


# consider skipping non hetsym ligands FIXME
def parseusercalculation(pdbname, session, debug=True, ignore_ligand_preset=False, ):
    logger = logging.getLogger('build')
    mypath = '/tmp/interactions/' + session + '/results/' + pdbname + '/output'
    module_dir = '/tmp/interactions/' + session
    results = []

    for f in listdir(mypath):
        if os.path.isfile(os.path.join(mypath, f)):
            result = yaml.load(open(mypath + "/" + f, 'rb'), Loader=yaml.FullLoader)
            output = result

            temp = f.replace('.yaml', '').split("_")
            temp.append([output])
            temp.append(round(output['score']))
            temp.append((output['inchikey']).strip())
            temp.append((output['smiles']).strip())
            results.append(temp)

            if 'prettyname' not in output:
                output['prettyname'] = temp[1]

    results = sorted(results, key=itemgetter(3), reverse=True)
    return results

# DEPRECATED
def showcalculation(request):

    context = calculate(request)

    return render(request, 'interaction/diagram.html', context)

# NOTE: this function is solely used by the sitesearch functionality
def calculate(request, redirect=None):
    if request.method == 'POST':
        form = PDBform(request.POST, request.FILES)
        if form.is_valid():

            pdbname = form.cleaned_data['pdbname'].strip()
            results = ''

            if not request.session.exists(request.session.session_key):
                request.session.create()
            session_key = request.session.session_key

            module_dirs = []
            module_dir = '/tmp/interactions'
            module_dirs.append(module_dir)
            module_dir = os.sep.join([module_dir, session_key])
            module_dirs.append(module_dir)
            module_dirs.append(os.sep.join([module_dir, 'pdbs']))
            module_dirs.append(os.sep.join([module_dir, 'temp']))

            # create dirs and set permissions (needed on some systems)
            for mdir in module_dirs:
                os.makedirs(mdir, exist_ok=True)
                os.chmod(mdir, 0o777)

            if 'file' in request.FILES:
                pdbdata = request.FILES['file']
                pdbname = os.path.splitext(str(pdbdata))[0]
                pdbname = pdbname.replace("_","")
                print(pdbname)
                with open(module_dir + '/pdbs/' + str(pdbdata).replace("_",""), 'wb+') as destination:
                    for chunk in pdbdata.chunks():
                        destination.write(chunk)

                temp_path = module_dir + '/pdbs/' + str(pdbdata).replace("_","")
                pdbdata = open(temp_path, 'r').read()
                runusercalculation(pdbname, session_key)

            else:
                pdbname = form.cleaned_data['pdbname'].strip()
                #print('pdbname selected!',pdbname)
                temp_path = module_dir + '/pdbs/' + pdbname + '.pdb'

                if not os.path.isfile(temp_path):
                    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbname
                    pdbdata = urllib.request.urlopen(
                        url).read().decode('utf-8')
                    f = open(temp_path, 'w')
                    f.write(pdbdata)
                    f.close()
                else:
                    pdbdata = open(temp_path, 'r').read()
                runusercalculation(pdbname, session_key)

            # MAPPING GPCRdb numbering onto pdb.
            generic_numbering = GenericNumbering(temp_path,top_results=1)
            out_struct = generic_numbering.assign_generic_numbers()
            structure_residues = generic_numbering.residues
            prot_id_list = generic_numbering.prot_id_list
            segments = {}

            generic_ids = []
            generic_number = []
            previous_seg = 'N-term'

            #Get segments built correctly for non-aligned residues
            for c, res in structure_residues.items():
                for i, r in sorted(res.items()):  # sort to be able to assign loops
                    if r.gpcrdb:
                        if r.gpcrdb[0] == '-':
                            # fix stefan - for bulge
                            r.gpcrdb = r.gpcrdb[1:] + "1"
                        r.gpcrdb = str(r.gpcrdb).replace('.', 'x')
                        generic_number.append(r.gpcrdb)
                    if r.gpcrdb_id:
                        generic_ids.append(r.gpcrdb_id)
                    if r.segment:
                        if not r.segment in segments:
                            segments[r.segment] = {}
                        segments[r.segment][r.number] = [
                            r.display, r.name, r.gpcrdb,r.residue_record]
                        previous_seg = r.segment
                    else:  # if no segment assigned by blast
                        if previous_seg in ['N-term', 'ICL1', 'ECL1', 'ICL2', 'ECL2', 'ICL3', 'ICL3', 'C-term']:
                            if not previous_seg in segments:
                                segments[previous_seg] = {}
                            segments[previous_seg][r.number] = ['', r.name, '',r.residue_record]
                        else:
                            if previous_seg == 'TM1':
                                previous_seg = 'ICL1'
                            elif previous_seg == 'TM2':
                                previous_seg = 'ECL1'
                            elif previous_seg == 'TM3':
                                previous_seg = 'ICL2'
                            elif previous_seg == 'TM4':
                                previous_seg = 'ECL2'
                            elif previous_seg == 'TM5':
                                previous_seg = 'ICL3'
                            elif previous_seg == 'TM6':
                                previous_seg = 'ECL3'
                            elif previous_seg == 'TM7':
                                previous_seg = 'C-term'
                            elif previous_seg == 'H8':
                                previous_seg = 'C-term'

                            if not previous_seg in segments:
                                segments[previous_seg] = {}
                            segments[previous_seg][r.number] = ['', r.name, '',r.residue_record]

            residue_list = []
            for seg, reslist in segments.items():
                for seq_number, v in sorted(reslist.items()):
                    if v[3]: #if blast assigned a residue, then use it. Otherwise just make something empty
                        r = v[3]
                    else:
                        r = Residue()
                    r.sequence_number = seq_number
                    r.segment_slug = seg
                    r.amino_acid = v[1]
                    residue_list.append(r)

            xtal = {}
            hetsyn = {}
            hetsyn_reverse = {}
            for line in pdbdata.splitlines():
                if line.startswith('HETSYN'):
                    # need to fix bad PDB formatting where col4 and col5 are
                    # put together for some reason -- usually seen when the id
                    # is +1000
                    m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)", line)
                    if (m):
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('HETNAM'):
                    # need to fix bad PDB formatting where col4 and col5 are
                    # put together for some reason -- usually seen when the id
                    # is +1000
                    m = re.match("HETNAM[\s]+([\w]{3})[\s]+(.+)", line)
                    if (m):
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('REVDAT   1'):
                    xtal['publication_date'] = line[13:22]
                    xtal['pdb_code'] = line[23:27]
                if line.startswith('JRNL        PMID'):
                    xtal['pubmed_id'] = line[19:].strip()
                if line.startswith('JRNL        DOI'):
                    xtal['doi_id'] = line[19:].strip()
                if line.startswith('REMARK   2 RESOLUTION.'):
                    xtal['resolution'] = line[22:].strip()

            results = parseusercalculation(
                pdbname, request.session.session_key)

            simple = collections.OrderedDict()
            simple_generic_number = collections.OrderedDict()
            residues_browser = []
            residue_table_list = []
            mainligand = ''

            for ligand in results:
                ligand_score = round(ligand[2][0]['score'])

                # select top hit
                if mainligand == '':
                    mainligand = ligand[1]

                simple[ligand[1]] = {'score': ligand_score}
                simple_generic_number[ligand[1]] = {'score': ligand_score}
                for interaction in ligand[2][0]['interactions']:
                    aa, pos, chain = regexaa(interaction[0])
                    if int(pos) in structure_residues[chain]:
                        r = structure_residues[chain][int(pos)]
                        display = r.display
                        segment = r.segment
                        generic = r.gpcrdb

                        if generic != "":
                            residue_table_list.append(generic)

                            if generic not in simple_generic_number[ligand[1]]:
                                simple_generic_number[ligand[1]][generic] = []
                            simple_generic_number[ligand[1]][generic].append(interaction[2])
                    else:
                        display = ''
                        segment = ''

                    if interaction[0] in simple[ligand[1]]:
                        simple[ligand[1]][interaction[0]].append(interaction[2])
                    else:
                        simple[ligand[1]][interaction[0]] = [interaction[2]]

                    residues_browser.append({'type': interaction[3], 'aa': aa, 'ligand': ligand[
                                            1], 'pos': pos, 'gpcrdb': display, 'segment': segment, 'slug':interaction[2]})
                break  # only use the top one

            # RESIDUE TABLE
            segments = ProteinSegment.objects.all().filter().prefetch_related()
            proteins = []
            protein_list = Protein.objects.filter(pk__in=prot_id_list)
            numbering_schemes_selection = [settings.DEFAULT_NUMBERING_SCHEME]
            for p in protein_list:
                proteins.append(p)
                if p.residue_numbering_scheme.slug not in numbering_schemes_selection:
                    numbering_schemes_selection.append(
                        p.residue_numbering_scheme.slug)

            numbering_schemes = ResidueNumberingScheme.objects.filter(
                slug__in=numbering_schemes_selection).all()
            default_scheme = numbering_schemes.get(
                slug=settings.DEFAULT_NUMBERING_SCHEME)
            data = OrderedDict()

            for segment in segments:
                data[segment.slug] = OrderedDict()
                residues = Residue.objects.filter(protein_segment=segment,  protein_conformation__protein__in=proteins,
                                                  generic_number__label__in=residue_table_list).prefetch_related('protein_conformation__protein',
                                                                                                                 'protein_conformation__state', 'protein_segment',
                                                                                                                 'generic_number', 'display_generic_number', 'generic_number__scheme',
                                                                                                                 'alternative_generic_numbers__scheme')
                for scheme in numbering_schemes:
                    if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                        for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                            data[segment.slug][pos] = {
                                scheme.slug: pos, 'seq': ['-'] * len(proteins)}
                    elif scheme == default_scheme:
                        for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                            data[segment.slug][pos] = {
                                scheme.slug: pos, 'seq': ['-'] * len(proteins)}

                for residue in residues:
                    alternatives = residue.alternative_generic_numbers.all()
                    pos = residue.generic_number
                    for alternative in alternatives:
                        if alternative.scheme not in numbering_schemes:
                            continue
                        scheme = alternative.scheme
                        if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                            pos = residue.generic_number
                            if scheme == pos.scheme:
                                data[segment.slug][pos.label]['seq'][proteins.index(
                                    residue.protein_conformation.protein)] = str(residue)
                            else:
                                if scheme.slug not in data[segment.slug][pos.label].keys():
                                    data[segment.slug][pos.label][
                                        scheme.slug] = alternative.label
                                if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                    data[segment.slug][pos.label][
                                        scheme.slug] += " " + alternative.label
                                data[segment.slug][pos.label]['seq'][proteins.index(
                                    residue.protein_conformation.protein)] = str(residue)
                        else:
                            if scheme.slug not in data[segment.slug][pos.label].keys():
                                data[segment.slug][pos.label][
                                    scheme.slug] = alternative.label
                            if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                data[segment.slug][pos.label][
                                    scheme.slug] += " " + alternative.label
                            data[segment.slug][pos.label]['seq'][proteins.index(
                                residue.protein_conformation.protein)] = str(residue)

            # Preparing the dictionary of list of lists. Dealing with tripple
            # nested dictionary in django templates is a nightmare
            flattened_data = OrderedDict.fromkeys(
                [x.slug for x in segments], [])
            for s in iter(flattened_data):
                flattened_data[s] = [[data[s][x][
                    y.slug] for y in numbering_schemes] + data[s][x]['seq'] for x in sorted(data[s])]

            context = {}
            context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [
                                    x.name for x in proteins], [x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
            context['segments'] = [
                x.slug for x in segments if len(data[x.slug])]
            context['data'] = flattened_data
            context['number_of_schemes'] = len(numbering_schemes)

            if redirect:
                # get simple selection from session
                simple_selection = request.session.get('selection', False)

                # create full selection and import simple selection (if it
                # exists)
                selection = Selection()
                if simple_selection:
                    selection.importer(simple_selection)

                # convert identified interactions to residue features and add them to the session
                # numbers in lists represent the interaction "hierarchy", i.e. if a residue has more than one
                # interaction,
                interaction_name_dict = {
                    'polar_double_neg_protein': [1, 'neg'],
                    'polar_double_pos_protein': [1, 'neg'],
                    'polar_pos_protein': [2, 'pos'],
                    'polar_neg_protein': [3, 'neg'],
                    'polar_neg_ligand': [4, 'hbd'],
                    'polar_pos_ligand': [5, 'hba'],
                    'polar_unknown_protein': [5, 'charge'],
                    'polar_donor_protein': [6, 'hbd'],
                    'polar_acceptor_protein': [7, 'hba'],
                    'polar_unspecified': [8, 'hb'],
                    'aro_ff': [9, 'ar'],
                    'aro_ef_protein': [10, 'ar'],
                    'aro_fe_protein':  [11, 'ar'],
                    'aro_ion_protein':  [12, 'pos'],
                    'aro_ion_ligand':  [12, 'ar'],
                }

                interaction_counter = 0
                for gn, interactions in simple_generic_number[mainligand].items():
                    if gn != 'score' and gn != 0.0:  # FIXME leave these out when dict is created
                        feature = False
                        for interaction in interactions:
                            if interaction in interaction_name_dict:
                                if (not feature
                                    or interaction_name_dict[interaction][0] < interaction_name_dict[feature][0]):
                                    feature = interaction

                        if not feature:
                            continue

                        # get residue number equivalent object
                        rne = ResidueGenericNumberEquivalent.objects.get(label=gn, scheme__slug='gpcrdba')

                        # create a selection item
                        properties = {
                            'feature': interaction_name_dict[feature][1],
                            'amino_acids': ','.join(definitions.AMINO_ACID_GROUPS_OLD[interaction_name_dict[feature][1]])
                        }
                        selection_item = SelectionItem(
                            'site_residue', rne, properties)

                        # add to selection
                        selection.add('segments', 'site_residue',
                                      selection_item)

                        # update the minimum match count for the active group
                        interaction_counter += 1
                        selection.site_residue_groups[selection.active_site_residue_group - 1][0] = interaction_counter

                # export simple selection that can be serialized
                simple_selection = selection.exporter()

                # add simple selection to session
                request.session['selection'] = simple_selection

                # re-direct to segment selection (with the extracted interactions already selected)
                return HttpResponseRedirect(redirect)
            else:
                # Only relevant when not redirecting - moved here
                HelixBox = DrawHelixBox(
                    residue_list, 'Class A', str('test'), nobuttons=1)
                SnakePlot = DrawSnakePlot(
                    residue_list, 'Class A', str('test'), nobuttons=1)

                return {'result': "Looking at " + pdbname, 'outputs': results,
                                                                    'simple': simple, 'simple_generic_number': simple_generic_number, 'xtal': xtal, 'pdbname': pdbname, 'mainligand': mainligand, 'residues': residues_browser,
                                                                    'HelixBox': HelixBox, 'SnakePlot': SnakePlot, 'data': context['data'],
                                                                    'header': context['header'], 'segments': context['segments'], 'number_of_schemes': len(numbering_schemes), 'proteins': proteins}

        else:
            print(form.errors)
            return HttpResponse("Error with form ")
    else:
        return HttpResponse("Ooops how did you get here?")


def download(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')

    session = request.GET.get('session')

    if session:
        session = request.session.session_key
        pdbdata = open('/tmp/interactions/' + session + '/results/' + pdbname +
                       '/interaction/' + pdbname + '_' + ligand + '.pdb', 'r').read()
        response = HttpResponse(pdbdata, content_type='text/plain')
    else:

        pair = StructureLigandInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(
            Q(ligand__properities__inchikey=ligand) | Q(ligand__name=ligand)).exclude(pdb_file__isnull=True).get()
        response = HttpResponse(pair.pdb_file.pdb, content_type='text/plain')
    return response

def excel(request, slug, **response_kwargs):
    if ('session' in response_kwargs):
        session = request.session.session_key

        mypath = '/tmp/interactions/' + session + '/pdbs/' + slug + '.pdb'

        generic_numbering = GenericNumbering(mypath)
        out_struct = generic_numbering.assign_generic_numbers()
        structure_residues = generic_numbering.residues
        results = parseusercalculation(slug,session)

        data = []
        for interaction in results[0][2][0]['interactions']:
            aa, pos, chain = regexaa(interaction[0])
            if int(pos) in structure_residues[chain]:
                r = structure_residues[chain][int(pos)]
                display = r.display
                segment = r.segment
                generic = r.gpcrdb
            else:
                display = ''
                segment = ''
            row = {}
            row['Sequence Number'] = pos
            row['Amino Acid'] = aa
            row['Generic Number'] = display
            row['Segment'] = segment
            row['Interaction'] = interaction[3]
            row['Interaction Slug'] = interaction[2]
            row['Ligand'] = results[0][2][0]['prettyname']

            data.append(row)

    else:

        interactions = ResidueFragmentInteraction.objects.filter(
            structure_ligand_pair__structure__pdb_code__index=slug, structure_ligand_pair__annotated=True).order_by('rotamer__residue__sequence_number')
        # return HttpResponse("Hello, world. You're at the polls index. "+slug)
        data = []
        for interaction in interactions:
            row = {}
            row['Sequence Number'] = interaction.rotamer.residue.sequence_number
            row['Amino Acid'] = interaction.rotamer.residue.amino_acid
            if interaction.rotamer.residue.display_generic_number:
                row['Generic Number'] = interaction.rotamer.residue.display_generic_number.label
                row['Segment'] = interaction.rotamer.residue.protein_segment.slug
            else:
                row['Generic Number'] = 'N/A'
                row['Segment'] = '-'

            row['Interaction'] = interaction.interaction_type.name
            row['Interaction Slug'] = interaction.interaction_type.slug
            row['Ligand'] = interaction.structure_ligand_pair.ligand.name

            data.append(row)


    headers = ['Ligand','Amino Acid','Sequence Number','Generic Number','Segment','Interaction','Interaction Slug']

    #EXCEL SOLUTION
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output)
    worksheet = workbook.add_worksheet()

    col = 0
    for h in headers:
        worksheet.write(0, col, h)
        col += 1
    row = 1
    for d in data:
        col = 0
        for h in headers:
            worksheet.write(row, col, d[h])
            col += 1
        row += 1
    workbook.close()
    output.seek(0)
    xlsx_data = output.read()

    response = HttpResponse(xlsx_data,content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Disposition'] = 'attachment; filename=Interaction_data_%s.xlsx' % slug
    return response

def ajax(request, slug, **response_kwargs):

    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment','display_generic_number','generic_number')
    lookup = {}

    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number.label] = r.sequence_number

    interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name=slug, structure_ligand_pair__annotated=True).exclude(interaction_type__type ='hidden').order_by('rotamer__residue__sequence_number')
    # return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for interaction in interactions:
        if interaction.rotamer.residue.generic_number:
            sequence_number = interaction.rotamer.residue.sequence_number
            if interaction.rotamer.residue.generic_number.label in lookup:
                sequence_number = lookup[interaction.rotamer.residue.generic_number.label]

            label = interaction.rotamer.residue.generic_number.label
            aa = interaction.rotamer.residue.amino_acid
            interactiontype = interaction.interaction_type.name
            if sequence_number not in jsondata:
                jsondata[sequence_number] = []
            jsondata[sequence_number].append([aa, interactiontype])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)


def ajaxLigand(request, slug, ligand, **response_kwargs):
    print(ligand)
    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment','display_generic_number','generic_number')
    lookup = {}

    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number.label] = r.sequence_number


    interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name=slug, structure_ligand_pair__ligand__name=ligand).exclude(interaction_type__type ='hidden').order_by('rotamer__residue__sequence_number')
    # return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for interaction in interactions:
        sequence_number = interaction.rotamer.residue.sequence_number
        sequence_number = lookup[interaction.rotamer.residue.generic_number.label]
        aa = interaction.rotamer.residue.amino_acid
        interactiontype = interaction.interaction_type.name
        if sequence_number not in jsondata:
            jsondata[sequence_number] = []
        jsondata[sequence_number].append([aa, interactiontype])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)


def pdbfragment(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    fragment = request.GET.get('fragment')

    result = ResidueFragmentInteraction.objects.filter(id=fragment).get()
    response = HttpResponse(result.rotamer.pdbdata.pdb +
                            result.fragment.pdbdata.pdb, content_type='text/plain')
    return response


def pdb(request):
    pdbname = request.GET.get('pdb')
    session = request.GET.get('session')
    # response = HttpResponse(mimetype='application/force-download')
    # #response['Content-Disposition'] = 'attachment; filename=%s' % smart_str(file_name)
    # response['Content-Disposition'] = 'attachment; filename=%s' % smart_str(pdbname+'_'+ligand+'.pdb')
    # mypath = module_dir+'/temp/results/'+pdbname+'/interaction/'+pdbname+'_'+ligand+'.pdb'
    # response['X-Sendfile'] = smart_str(mypath)
    if session:
        session = request.session.session_key
        pdbdata = open('/tmp/interactions/' + session +
                       '/pdbs/' + pdbname + '.pdb', 'r').read()
        response = HttpResponse(pdbdata, content_type='text/plain')
    else:
        web_resource = WebResource.objects.get(slug='pdb')
        web_link, created = WebLink.objects.get_or_create(
            web_resource=web_resource, index=pdbname)

        structure = Structure.objects.filter(pdb_code=web_link)
        if structure.exists():
            structure = Structure.objects.get(pdb_code=web_link)
        else:
            quit()  # quit!

        if structure.pdb_data is None:
            quit()

        response = HttpResponse(structure.pdb_data.pdb,
                                content_type='text/plain')
    return response
