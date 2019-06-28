from django.conf import settings
from django.shortcuts import render
from django.db.models import Count, Avg, Min, Max, Q
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView, View

from interaction.models import ResidueFragmentInteraction
from mutation.models import MutationExperiment
from protein.models import Protein
from residue.models import Residue
from structure.models import Structure

def hotspotsView(request):
    """
    Show hotspots viewer page
    """
    return render(request, 'hotspots/hotspotsView.html')


def getHotspots(request):
    data = {'error': 0}

    # DEBUG: for now only Class A
    gpcr_class = "001"

    # Obtain all proteins
    class_proteins = Protein.objects.filter(family__slug__startswith=gpcr_class, sequence_type__slug='wt', species__common_name='Human')
    class_count = class_proteins.count()

    # SEQUENCE: number of same amino acids per position in class
    residues = Residue.objects.filter(protein_conformation__protein__in = class_proteins)\
                    .exclude(generic_number=None)\
                    .values('generic_number__label','amino_acid')\
                    .annotate(number_occurrences=Count("protein_conformation"))\
                    .order_by('generic_number__label')

    seq_conservation = {}
    for entry in list(residues):
        if not entry["generic_number__label"] in seq_conservation:
            seq_conservation[entry["generic_number__label"]] = {}
        seq_conservation[entry["generic_number__label"]][entry["amino_acid"]] = entry["number_occurrences"]

    data["seq_conservation"] = seq_conservation

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

    data["mutation_count"] = mutation_count

    # STRUCTURE: Ligand contacts per position per protein (consider also per subfamily/ligand type/class)
    ligand_interactions = ResidueFragmentInteraction.objects.filter(\
                structure_ligand_pair__structure__protein_conformation__protein__parent__in=class_proteins)\
                .exclude(interaction_type__type='hidden')\
                .exclude(rotamer__residue__generic_number=None)\
                .values("rotamer__residue__protein_conformation__protein__parent__entry_name", "rotamer__residue__generic_number__label")\
                .annotate(unique_contacts=Count("rotamer__structure", distinct=True))\
                .order_by('rotamer__residue__protein_conformation__protein__parent__entry_name', "rotamer__residue__generic_number__label")

    contact_count = {}
    for entry in list(ligand_interactions):
        if not entry["rotamer__residue__protein_conformation__protein__parent__entry_name"] in contact_count:
            contact_count[entry["rotamer__residue__protein_conformation__protein__parent__entry_name"]] = {}
        contact_count[entry["rotamer__residue__protein_conformation__protein__parent__entry_name"]][entry["rotamer__residue__generic_number__label"]] = entry["unique_contacts"]

    data["contact_count"] = contact_count

    return JsonResponse(data)
