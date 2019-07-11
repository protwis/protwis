from django.conf import settings
from django.shortcuts import render
from django.db.models import Count, Avg, Min, Max, Q
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView, View

from interaction.models import ResidueFragmentInteraction
from mutation.models import MutationExperiment
from protein.models import Protein, ProteinSegment
from residue.models import Residue
from structure.models import Structure

from collections import OrderedDict
import functools

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

def hotspotsView(request):
    """
    Show hotspots viewer page
    """
    return render(request, 'hotspots/hotspotsView.html')

# @cache_page(60*60*24*7)
def getHotspots(request):
    def gpcrdb_number_comparator(e1, e2):
        t1 = e1.split('x')
        t2 = e2.split('x')

        if e1 == e2:
            return 0

        if t1[0] == t2[0]:
            if t1[1] < t2[1]:
                return -1
            else:
                return 1

        if t1[0] < t2[0]:
            return -1
        else:
            return 1

    data = {'error': 0}

    # DEBUG: for now only Class A
    gpcr_class = "001"

    # Obtain all proteins
    class_proteins = Protein.objects.filter(family__slug__startswith=gpcr_class, sequence_type__slug='wt', species__common_name='Human')\
                            .prefetch_related('family__parent','family__parent__parent')
    class_count = class_proteins.count()

    protein_dictionary = {}
    for p in class_proteins:
        protein_dictionary[p.entry_name] = {}
        protein_dictionary[p.entry_name]["receptor_family"] = p.family.parent.short()
        protein_dictionary[p.entry_name]["ligand_type"] = p.family.parent.parent.short()

    # SEQUENCE: number of same amino acids per position in class
    residues = Residue.objects.filter(protein_conformation__protein__in = class_proteins)\
                    .exclude(generic_number=None)\
                    .values('generic_number__label','amino_acid')\
                    .annotate(number_occurrences=Count("protein_conformation"))\
                    .order_by('generic_number__label') # Necessary otherwise the key is used -> messing up the count

    seq_conservation = {}
    for entry in list(residues):
        if not entry["generic_number__label"] in seq_conservation:
            seq_conservation[entry["generic_number__label"]] = {}
        seq_conservation[entry["generic_number__label"]][entry["amino_acid"]] = entry["number_occurrences"]

    #data["seq_conservation"] = seq_conservation

    # MUTATION: obtain ligand mutations >5 fold effect
    receptor_mutations = MutationExperiment.objects.filter(Q(foldchange__gte = 5) | Q(foldchange__lte = -5), protein__in = class_proteins)\
                    .exclude(residue__generic_number=None)\
                    .values("protein__entry_name", "residue__generic_number__label")\
                    .annotate(unique_mutations=Count("id"))\
                    .order_by('protein__entry_name', "residue__generic_number__label")

    mutation_count = {}
    for entry in list(receptor_mutations):
        if not entry["protein__entry_name"] in mutation_count:
            mutation_count[entry["protein__entry_name"]] = {}
        mutation_count[entry["protein__entry_name"]][entry["residue__generic_number__label"]] = entry["unique_mutations"]

    #data["mutation_count"] = mutation_count

    # STRUCTURE: Ligand contacts per position per protein (consider also per subfamily/ligand type/class)
    ligand_interactions = ResidueFragmentInteraction.objects.filter(\
                structure_ligand_pair__structure__protein_conformation__protein__parent__in=class_proteins)\
                .exclude(structure_ligand_pair__annotated=False)\
                .exclude(rotamer__residue__generic_number=None)\
                .exclude(rotamer__residue__generic_number=None)\
                .values("rotamer__residue__protein_conformation__protein__parent__entry_name", "rotamer__residue__generic_number__label")\
                .annotate(unique_contacts=Count("rotamer__structure", distinct=True))\
                .order_by('rotamer__residue__protein_conformation__protein__parent__entry_name', "rotamer__residue__generic_number__label")

    contact_count = {}
    for entry in list(ligand_interactions):
        if not entry["rotamer__residue__protein_conformation__protein__parent__entry_name"] in contact_count:
            contact_count[entry["rotamer__residue__protein_conformation__protein__parent__entry_name"]] = {}
        contact_count[entry["rotamer__residue__protein_conformation__protein__parent__entry_name"]][entry["rotamer__residue__generic_number__label"]] = entry["unique_contacts"]

    #data["contact_count"] = contact_count

    # MISSING: alignment per entry
    aln = Alignment()
    aln.load_proteins(class_proteins)

    # refine later
    aln.load_segments(ProteinSegment.objects.filter(slug__in=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']))
    aln.build_alignment()

    # start parsing the data
    residue_matrix = {}
    generic_numbers = set()
    for i, p in enumerate(aln.unique_proteins):
        entry_name = p.protein.entry_name

        residue_matrix[entry_name] = protein_dictionary[entry_name]
        for j, s in p.alignment.items():
            for p in s:
                generic_number = p[0]
                generic_numbers.add(generic_number)
                display_generic_number = p[1]
                amino_acid = p[2]

                # TODO: handle gaps

                # OTPIMIZE THIS by zipping efficiently
                seq_count = 0
                if generic_number in seq_conservation and amino_acid in seq_conservation[generic_number]:
                    seq_count = seq_conservation[generic_number][amino_acid]

                mut_count = 0
                if entry_name in mutation_count and generic_number in mutation_count[entry_name]:
                    mut_count = mutation_count[entry_name][generic_number]

                con_count = 0
                if entry_name in contact_count and generic_number in contact_count[entry_name]:
                    con_count = contact_count[entry_name][generic_number]

                residue_matrix[entry_name][generic_number] = [amino_acid, display_generic_number, ['#f0fcfa',seq_count], ['#fbf0fc',mut_count], ['#ccc',con_count]]

    data['sorted_gns'] = sorted(generic_numbers, key=functools.cmp_to_key(gpcrdb_number_comparator))
    data["residue_matrix"] = residue_matrix

    return JsonResponse(data)
