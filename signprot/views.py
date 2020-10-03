from django.contrib.postgres.aggregates import ArrayAgg
from django.core.cache import cache
from django.db.models import F, Q
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView

from common import definitions
from common.diagrams_gpcr import DrawSnakePlot
from common.diagrams_gprotein import DrawGproteinPlot
from common.tools import fetch_from_web_api
from common.views import AbsTargetSelection
from contactnetwork.models import InteractingResiduePair
from mutation.models import MutationExperiment
from protein.models import (Gene, Protein, ProteinAlias, ProteinConformation, ProteinFamily, ProteinGProtein,
                            ProteinGProteinPair, ProteinSegment)
from residue.models import (Residue, ResidueGenericNumberEquivalent, ResiduePositionSet)
from seqsign.sequence_signature import (SequenceSignature, SignatureMatch)
from signprot.interactions import (get_entry_names, get_generic_numbers, get_ignore_info, get_protein_segments,
                                   get_signature_features, group_signature_features, prepare_signature_match)
from signprot.models import (SignprotBarcode, SignprotComplex, SignprotStructure)
from structure.models import Structure

import json
import re
import time
from collections import Counter, OrderedDict
from decimal import Decimal
from pprint import pprint


class BrowseSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    psets = False
    filters = True
    filter_gprotein = True

    type_of_selection = 'browse_gprot'

    description = 'Select a G protein or family by searching or browsing in the right column.'
    description = 'Select a G protein (family) by searching or browsing in the middle. The selection is viewed to' \
                  + ' the right.'
    docs = 'receptors.html'
    target_input = False

    selection_boxes = OrderedDict([
        ('reference', False), ('targets', True),
        ('segments', False),
    ])
    try:
        ppf_g = ProteinFamily.objects.get(slug="100_001")
        # ppf_a = ProteinFamily.objects.get(slug="200_000")
        # pfs = ProteinFamily.objects.filter(parent__in=[ppf_g.id,ppf_a.id])
        pfs = ProteinFamily.objects.filter(parent__in=[ppf_g.id])
        ps = Protein.objects.filter(family__in=[ppf_g])  # ,ppf_a
        tree_indent_level = []
        # action = 'expand'
        # remove the parent family (for all other families than the root of the tree, the parent should be shown)
        # del ppf_g
        # del ppf_a
    except Exception as e:
        pass


class ArrestinSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    psets = False
    filters = True
    filter_gprotein = True

    type_of_selection = 'browse_gprot'

    description = 'Select an Arrestin (family) by searching or browsing in the middle. The selection is viewed to' \
                  + ' the right.'
    docs = 'signalproteins.html'
    target_input = False

    selection_boxes = OrderedDict([
        ('reference', False), ('targets', True),
        ('segments', False),
    ])
    try:
        if ProteinFamily.objects.filter(slug="200_000").exists():
            ppf = ProteinFamily.objects.get(slug="200_000")
            pfs = ProteinFamily.objects.filter(parent=ppf.id)
            ps = Protein.objects.filter(family=ppf)

            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass


class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    filters = False
    psets = False
    target_input = False
    redirect_on_select = True
    type_of_selection = 'ginterface'
    title = 'SELECT TARGET for Gs INTERFACE'
    description = 'Select a reference target by searching or browsing.' \
                  + '\n\nThe Gs interface from adrb2 (PDB: 3SN6) will be superposed onto the selected target.' \
                  + '\n\nAn interaction browser for the adrb2 Gs interface will be given for comparison"'

    # template_name = 'common/targetselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '#',
            'color': 'success',
        },
    }


class CouplingBrowser(TemplateView):
    """
    Class based generic view which serves coupling data between Receptors and G-proteins.
    Data coming from Guide to Pharmacology, Asuka Inuoue and Michel Bouvier.
    More data might come later from Roth and Strachan TRUPATH biosensor.
    :param dataset: ProteinGProteinPair (see build/management/commands/build_g_proteins.py)
    :return: context
    """
    template_name = "signprot/coupling_browser.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        fam_tab = self.fams_tab()

        protobj = Protein.objects.filter(sequence_type__slug='wt',
                                         family__slug__startswith='00',
                                         species__common_name='Human').prefetch_related(
            "family",
            "parent",
            "source")

        coupobj = ProteinGProteinPair.objects.all().prefetch_related('protein',
                                                                     'g_protein_subunit',
                                                                     'g_protein')

#        context['proteins'] = protobj
#        context['couplings'] = coupobj
        context['famtab'] = fam_tab

        return context

    @staticmethod
    def fams_tab():
        """
        This function returns what is needed in the families (maybe others - needs refactoring) tab of the coupling
        browser.
        :return: key.value pairs from fd dictionary
        keys = gene names of proteins
        values = gpcrDB general values and emax-double-normalized, pec50-double-normalized, emax-mean, pec50-mean
        emax-sem, pec50sem.
        """
        threshold_primary_inoue = 0.5  # -0.1
        threshold_secondary_inoue = 0.01  # -1
        threshold_primary_bouvier = 0.5  # -0.1
        threshold_secondary_bouvier = 0.01  # -1
        proteins = Protein.objects.filter(sequence_type__slug='wt',
                                          family__slug__startswith='00',
                                          species__common_name='Human').all().prefetch_related('family')
        data = {}
        class_names = {}
        family_names = {}
        for p in proteins:
            p_class = p.family.slug.split('_')[0]
            if p_class not in class_names:
                class_names[p_class] = p.family.parent.parent.parent.name
                family_names[p_class] = p.family.parent.name
            p_class_name = class_names[p_class].replace('Class', '').strip()
            p_family_name = family_names[p_class].replace('receptors', '').strip()
            p_accession = p.accession
            p_entryname = p.entry_name
            data[p.entry_short()] = {'class': p_class_name,
                                     'family': p_family_name,
                                     'accession': p_accession,
                                     'entryname': p_entryname,
                                     'pretty': p.short(),
                                     'GuideToPharma': {},
                                     'Aska': {},
                                     'Bouvier': {}}

        distinct_g_families = []
        distinct_g_subunit_families = {}
        couplings = ProteinGProteinPair.objects.prefetch_related('protein',
                                                                 'g_protein_subunit',
                                                                 'g_protein')
        for c in couplings:
            p = c.protein.entry_short()
            
            # Skip entries without any annotation
            if p not in data:
                continue

            s = c.source
            t = c.transduction
            emdn = c.emax_dnorm
            pec50dn = c.pec50_dnorm
            emmean = c.emax_mean
            pec50mean = c.pec50_mean
            emsem = c.emax_sem
            pec50sem = c.pec50_sem
            gf = c.g_protein.name
            gf = gf.replace(" family", "")

            if gf not in distinct_g_families:
                distinct_g_families.append(gf)
                distinct_g_subunit_families[gf] = []

            if c.g_protein_subunit:
                g = c.g_protein_subunit.entry_name
                g = g.replace("_human", "")
                if g not in distinct_g_subunit_families[gf]:
                    distinct_g_subunit_families[gf].append(g)
                    distinct_g_subunit_families[gf] = sorted(distinct_g_subunit_families[gf])

            if s not in data[p]:
                data[p][s] = {}

            if gf not in data[p][s]:
                data[p][s][gf] = {}

            # If transduction (primary, secondary) in Guide To Pharmacology data
            if t:
                data[p][s][gf] = t

            # Else t (primary secondary classification) doesn't exist, such as
            # source [s]  Aska (Inoue) or Bouvier.
            else:
                if 'emdn' not in data[p][s][gf]:
                    data[p][s][gf] = {'emdn': {}, 'max': 0.00, 'maxid': str()}
                if emdn is None:
                    data[p][s][gf]['emdn'][g] = round(Decimal(0), 2)
                    continue
                data[p][s][gf]['emdn'][g] = round(Decimal(emdn), 2)
                # get the highest number into 'max'
                # TODO: Who is max? It will eventually be needed, count on it!
                if emdn > data[p][s][gf]['max']:
                    data[p][s][gf]['max'] = round(Decimal(emdn), 2)
                    data[p][s][gf]['maxid'] = g
                # pec50dn
                if 'pec50dn' not in data[p][s][gf]:
                    data[p][s][gf]['pec50dn'] = {}
                if pec50dn is None:
                    data[p][s][gf]['pec50dn'][g] = round(Decimal(0), 2)
                    continue
                data[p][s][gf]['pec50dn'][g] = round(Decimal(pec50dn), 2)
                # emmean
                if 'emmean' not in data[p][s][gf]:
                    data[p][s][gf]['emmean'] = {}
                if emmean is None:
                    data[p][s][gf]['emmean'][g] = round(Decimal(0), 2)
                    continue
                data[p][s][gf]['emmean'][g] = round(Decimal(emmean), 2)
                # pec50mean
                if 'pec50mean' not in data[p][s][gf]:
                    data[p][s][gf]['pec50mean'] = {}
                if pec50mean is None:
                    data[p][s][gf]['pec50mean'][g] = round(Decimal(0), 2)
                    continue
                data[p][s][gf]['pec50mean'][g] = round(Decimal(pec50mean), 2)
                # emsem
                if 'emsem' not in data[p][s][gf]:
                    data[p][s][gf]['emsem'] = {}
                if emsem is None:
                    data[p][s][gf]['emsem'][g] = round(Decimal(0), 2)
                    continue
                data[p][s][gf]['emsem'][g] = round(Decimal(emsem), 2)
                # pec50sem
                if 'pec50sem' not in data[p][s][gf]:
                    data[p][s][gf]['pec50sem'] = {}
                if pec50sem is None:
                    data[p][s][gf]['pec50sem'][g] = round(Decimal(0), 2)
                    continue
                data[p][s][gf]['pec50sem'][g] = round(Decimal(pec50sem), 2)


        fd = {}  # final data
        distinct_g_families = ['Gs', 'Gi/Go', 'Gq/G11', 'G12/G13']
        distinct_g_subunit_families = OrderedDict(
            [('Gs', ['gnas2', 'gnal']),
             ('Gi/Go', ['gnai1', 'gnai2', 'gnai3', 'gnat1', 'gnat2', 'gnat3', 'gnao', 'gnaz']),
             ('Gq/G11', ['gnaq', 'gna11', 'gna14', 'gna15']),
             ('G12/G13', ['gna12', 'gna13'])])

        # This for loop, does the job of putting together
        # data-sets, GuideToPharma, Asuka's and Bouvier.
        for p, v in data.items():
            fd[p] = [v['class'], v['family'], v['accession'], v['entryname'], p, v['pretty']]
            s = 'GuideToPharma'
            # MERGED CRITERIA FOR COUPLING
            # merge primary, secondary, coupling, no coupling classificationm currently from
            # three sources, GtP, Inoue, Bouvier. Exact rules being worked on.
            for gf in distinct_g_families:
                values = []
                if 'GuideToPharma' in v and gf in v['GuideToPharma']:
                    values.append(v['GuideToPharma'][gf])
                else:
                    values.append('na gtf')
                if 'Aska' in v and gf in v['Aska']:
                    max = v['Aska'][gf]['max']
                    if max > threshold_primary_inoue:
                        values.append('primary')
                    elif max > threshold_secondary_inoue:
                        values.append('secondary')
                    elif max == 0:
                        values.append('NCI')
                else:
                    values.append('na inoue')
                if 'Bouvier' in v and gf in v['Bouvier']:
                    max = v['Bouvier'][gf]['max']
                    if max > threshold_primary_bouvier:
                        values.append('primary')
                    elif max > threshold_secondary_bouvier:
                        values.append('secondary')
                    elif max == 0:
                        values.append('NCB')
                else:
                    values.append('na bouvier')
                if 'primary' in values:
                    fd[p].append('primary')
                elif 'secondary' in values:
                    fd[p].append('secondary')
#                elif 'NCI' or 'NCB' in values:
#                    fd[p].append('no coupling')
                else:
                    fd[p].append('NA')  # Should this be NA of no-coupling? What's the GtP meaning?
                # if (len(values) >= 2):
                #     print(values, gf, p)

            # Loop over GuideToPharma
            s = 'GuideToPharma'
            for gf in distinct_g_families:
                if gf in v[s]:
                    fd[p].append(v[s][gf])
                else:
                    #fd[p].append("NAg")
                    fd[p].append('')

            # Loop over Aska
            s = 'Aska'
            for gf in distinct_g_families:
                if gf in v[s]:
                    if v[s][gf]['max'] > threshold_primary_inoue:
                        fd[p].append("primaryino")
                    elif v[s][gf]['max'] > threshold_secondary_inoue:
                        fd[p].append("secondaryino")
                    else:
                        fd[p].append("No coup. Ino")
                else:
                    fd[p].append('NAino')  # This last empty one means not-available, NOT no-coupling

            # Last bit adds values in subfamilies (i.e. subtypes)
            for gf, sfs in distinct_g_subunit_families.items():
                for sf in sfs:
                    if gf in v[s]:
                        if sf in v[s][gf]['emdn']:
                            fd[p].append(v[s][gf]['emdn'][sf])
                        else:
                            #fd[p].append("NAi1")
                            fd[p].append('')
                    else:
                        #fd[p].append("emdn nai")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['pec50dn']:
                            fd[p].append(v[s][gf]['pec50dn'][sf])
                        else:
                            #fd[p].append("NAi1")
                            fd[p].append('')
                    else:
                        #fd[p].append("pec50dn nai")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['emmean']:
                            fd[p].append(v[s][gf]['emmean'][sf])
                        else:
                            #fd[p].append("NAi1")
                            fd[p].append('')
                    else:
                        #fd[p].append("emmean nai")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['pec50mean']:
                            fd[p].append(v[s][gf]['pec50mean'][sf])
                        else:
                            #fd[p].append("NAi1")
                            fd[p].append('')
                    else:
                        #fd[p].append("pec50mean nai")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['emsem']:
                            fd[p].append(v[s][gf]['emsem'][sf])
                        else:
                           #fd[p].append("NAi1")
                           fd[p].append('')
                    else:
                        #fd[p].append("emsem nai")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['pec50sem']:
                            fd[p].append(v[s][gf]['pec50sem'][sf])
                        else:
                            #fd[p].append("NAi1")
                            fd[p].append('')
                    else:
                        #fd[p].append("pec50sem nai")
                        fd[p].append('')

            # Loop over Bouvier
            s = 'Bouvier'
            for gf in distinct_g_families:
                if gf in v[s]:
                    if v[s][gf]['max'] > threshold_primary_bouvier:
                        fd[p].append("primary")
                    elif v[s][gf]['max'] > threshold_secondary_bouvier:
                        fd[p].append("secondary")
                    else:
                        fd[p].append("NC")
                else:
                    fd[p].append('NA')  # This last empty one means not-available, NOT no-coupling

            # Last bit adds values to subfamilies (i.e. subtypes)
            for gf, sfs in distinct_g_subunit_families.items():
                for sf in sfs:
                    if gf in v[s]:
                        if sf in v[s][gf]['emdn']:
                            fd[p].append(v[s][gf]['emdn'][sf])
                        else:
                            #fd[p].append("NAb1")
                            fd[p].append('')
                    else:
                        #fd[p].append("emdn nab")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['pec50dn']:
                            fd[p].append(v[s][gf]['pec50dn'][sf])
                        else:
                            #fd[p].append("NAb1")
                            fd[p].append('')
                    else:
                        #fd[p].append("pec50dn nab")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['emmean']:
                            fd[p].append(v[s][gf]['emmean'][sf])
                        else:
                            #fd[p].append("NAb1")
                            fd[p].append('')
                    else:
                        #fd[p].append("emmean nab")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['pec50mean']:
                            fd[p].append(v[s][gf]['pec50mean'][sf])
                        else:
                            #fd[p].append("NAb1")
                            fd[p].append('')
                    else:
                        #fd[p].append("pec50mean nab")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['emsem']:
                            fd[p].append(v[s][gf]['emsem'][sf])
                        else:
                            #fd[p].append("NAb1")
                            fd[p].append('')
                    else:
                        #fd[p].append("emsem nab")
                        fd[p].append('')

                    if gf in v[s]:
                        if sf in v[s][gf]['pec50sem']:
                            fd[p].append(v[s][gf]['pec50sem'][sf])
                        else:
                            #fd[p].append("NAb1")
                            fd[p].append('')
                    else:
                        #fd[p].append("pec50sem nab")
                        fd[p].append('')

            # max Values for Gs, Gi/Go, Gq/G11, G12/13 and source Inoue
            s = 'Aska'
            for gf in distinct_g_families:
                if gf in v[s]:
                    # fd[p].append(v[s][gf]['max'] + 100)
                    fd[p].append(v[s][gf]['max'])
                else:
                    #fd[p].append("NAi max")
                    fd[p].append('')

            # max Values for Gs, Gi/Go, Gq/G11, G12/13 and source Bouvier
            s = 'Bouvier'
            for gf in distinct_g_families:
                if gf in v[s]:
                    #fd[p].append(v[s][gf]['max'] + 200)
                    fd[p].append(v[s][gf]['max'])
                else:
                    #fd[p].append("NAb max")
                    fd[p].append('')

        return fd

def GProtein(request, dataset="GuideToPharma"):
    name_of_cache = 'gprotein_statistics_{}'.format(dataset)

    context = cache.get(name_of_cache)
    if context == None:

        context = OrderedDict()
        i = 0
        gproteins = ProteinGProtein.objects.all().prefetch_related('proteingproteinpair_set')

        slug_translate = {'001': "ClassA", '002': "ClassB1", '004': "ClassC", '006': "ClassF"}
        selectivitydata = {}
        for slug in slug_translate.keys():
            jsondata = {}
            for gp in gproteins:
                # ps = gp.proteingproteinpair_set.all()
                ps = gp.proteingproteinpair_set.filter(protein__family__slug__startswith=slug,
                                                       source=dataset).prefetch_related('protein')
                # print(ps,len(ps))
                if ps:
                    jsondata[str(gp)] = []
                    for p in ps:
                        if dataset == "Aska" and p.log_rai_mean < -1:
                            continue
                        if str(p.protein.entry_name).split('_')[0].upper() not in selectivitydata:
                            selectivitydata[str(p.protein.entry_name).split('_')[0].upper()] = []
                        selectivitydata[str(p.protein.entry_name).split('_')[0].upper()].append(str(gp))
                        # print(p.protein.family.parent.parent.parent)
                        jsondata[str(gp)].append(str(p.protein.entry_name) + '\n')

                    jsondata[str(gp)] = ''.join(jsondata[str(gp)])

            context[slug_translate[slug]] = jsondata

        context["selectivitydata"] = selectivitydata

    cache.set(name_of_cache, context, 60 * 60 * 24 * 7)  # seven days timeout on cache

    return render(request,
                  'signprot/gprotein.html',
                  context
    )

# @cache_page(60*60*24*2) # 2 days caching
def couplings(request, template_name='signprot/coupling_browser.html'):
    """
    Presents coupling data between Receptors and G-proteins.
    Data coming from Guide to Pharmacology, Asuka Inuoue and Michel Bouvier
    """
    context = OrderedDict()
    threshold_primary = 0.5 # -0.1
    threshold_secondary = 0.01 # -1
    proteins = Protein.objects.filter(sequence_type__slug='wt',
                                      family__slug__startswith='00',
                                      species__common_name='Human').all().prefetch_related('family')
    data = {}
    class_names = {}
    family_names = {}

    for p in proteins:
        p_class = p.family.slug.split('_')[0]
        if p_class not in class_names:
            class_names[p_class] = p.family.parent.parent.parent.name
            family_names[p_class] = p.family.parent.name
        p_class_name = class_names[p_class].replace('Class','').strip()
        p_family_name = family_names[p_class].replace('receptors','').strip()
        p_accession = p.accession
        data[p.entry_short()] = {'class': p_class_name,
                                 'family': p_family_name,
                                 'accession': p_accession,
                                 'pretty': p.short(),
                                 'GuideToPharma': {},
                                 'Aska': {},
                                 'Bouvier': {}}
    distinct_g_families = []
    distinct_g_subunit_families = {}
    distinct_sources = ['GuideToPharma', 'Aska', 'Bouvier']
    couplings = ProteinGProteinPair.objects.all().prefetch_related('protein',
                                                                   'g_protein_subunit',
                                                                   'g_protein')

    for c in couplings:
        p = c.protein.entry_short()
        s = c.source
        t = c.transduction
        m = c.emax_dnorm
        gf = c.g_protein.name
        gf = gf.replace(" family", "")

        if gf not in distinct_g_families:
            distinct_g_families.append(gf)
            distinct_g_subunit_families[gf] = []

        if c.g_protein_subunit:
            g = c.g_protein_subunit.entry_name
            g = g.replace("_human", "")
            if g not in distinct_g_subunit_families[gf]:
                distinct_g_subunit_families[gf].append(g)
                distinct_g_subunit_families[gf] = sorted(distinct_g_subunit_families[gf])

        if s not in data[p]:
            data[p][s] = {}

        if gf not in data[p][s]:
            data[p][s][gf] = {}

        # If transduction in GuideToPharma data
        if t:
            data[p][s][gf] = t
        else:
            if 'subunits' not in data[p][s][gf]:
                data[p][s][gf] = {'subunits': {}, 'best': 0.00}
            if m is None:
                continue
            data[p][s][gf]['subunits'][g] = round(Decimal(m), 2)
            if round(Decimal(m), 2) == -0.00:
                data[p][s][gf]['subunits'][g] = 0.00
            # get the lowest number into 'best'
            if m > data[p][s][gf]['best']:
                data[p][s][gf]['best'] = round(Decimal(m), 2)
    fd = {}  # final data
#    distinct_g_families = sorted(distinct_g_families)
    distinct_g_families = ['Gs', 'Gi/Go', 'Gq/G11', 'G12/G13']
    distinct_g_subunit_families = OrderedDict(
        [('Gs', ['gnas2', 'gnal']),
         ('Gi/Go', ['gnai1', 'gnai2', 'gnai3', 'gnao', 'gnaz']),
         ('Gq/G11', ['gnaq', 'gna11', 'gna14', 'gna15']),
         ('G12/G13', ['gna12', 'gna13'])])

    # This for loop, which perhaps should be a function in itself, perhaps an instance of a Couplings class, does
    # the job of merging together two data-sets, that of the GuideToPharma and Asuka's results.
    for p, v in data.items():
        fd[p] = [v['class'], v['family'], v['accession'], p, v['pretty']]
        s = 'GuideToPharma'
        # Merge
        for gf in distinct_g_families:
            values = []
            if 'GuideToPharma' in v and gf in v['GuideToPharma']:
                values.append(v['GuideToPharma'][gf])
            if 'Aska' in v and gf in v['Aska']:
                best = v['Aska'][gf]['best']
                if best > threshold_primary:
                    values.append('primary')
                elif best > threshold_secondary:
                    values.append('secondary')
            if 'Bouvier' in v and gf in v['Bouvier']:
                best = v['Bouvier'][gf]['best']
                if best > threshold_primary:
                    values.append('primary')
                elif best > threshold_secondary:
                    values.append('secondary')
            if 'primary' in values:
                fd[p].append('primary')
            elif 'secondary' in values:
                fd[p].append('secondary')
            else:
                fd[p].append('')
        s = 'GuideToPharma'
        # First loop over GuideToPharma
        for gf in distinct_g_families:
            if gf in v[s]:
                fd[p].append(v[s][gf])
            else:
                fd[p].append("")
        s = 'Aska'
        for gf in distinct_g_families:
            if gf in v[s]:
                if v[s][gf]['best'] > threshold_primary:
                    fd[p].append("primary")
                elif v[s][gf]['best'] > threshold_secondary:
                    fd[p].append("secondary")
                else:
                    fd[p].append("No coupling")
            else:
                fd[p].append("")
        s = 'Bouvier'
        for gf in distinct_g_families:
            if gf in v[s]:
                if v[s][gf]['best'] > threshold_primary:
                    fd[p].append("primary")
                elif v[s][gf]['best'] > threshold_secondary:
                    fd[p].append("secondary")
                else:
                    fd[p].append("No coupling")
            else:
                fd[p].append("")
        for gf, sfs in distinct_g_subunit_families.items():
            for sf in sfs:
                if gf in v[s]:
                    if sf in v[s][gf]['subunits']:
                        fd[p].append(v[s][gf]['subunits'][sf])
                    else:
                        fd[p].append("")
                else:
                    fd[p].append("")

    context['data'] = fd
    context['distinct_gf'] = distinct_g_families
    context['distinct_sf'] = distinct_g_subunit_families

    return render(request,
                  template_name, context
    )

@cache_page(60 * 60 * 24 * 2)
def familyDetail(request, slug):
    # get family
    pf = ProteinFamily.objects.get(slug=slug)
    # get family list
    ppf = pf
    families = [ppf.name]
    while ppf.parent.parent:
        families.append(ppf.parent.name)
        ppf = ppf.parent
    families.reverse()

    # number of proteins
    proteins = Protein.objects.filter(family__slug__startswith=pf.slug, sequence_type__slug='wt')
    no_of_proteins = proteins.count()
    no_of_human_proteins = Protein.objects.filter(family__slug__startswith=pf.slug, species__id=1,
                                                  sequence_type__slug='wt').count()
    list_proteins = list(proteins.values_list('pk', flat=True))

    # get structures of this family
    structures = SignprotStructure.objects.filter(protein__family__slug__startswith=slug)
    complex_structures = SignprotComplex.objects.filter(protein__family__slug__startswith=slug)

    mutations = MutationExperiment.objects.filter(protein__in=proteins).prefetch_related('residue__generic_number',
                                                                                         'exp_qual', 'ligand')

    mutations_list = {}
    for mutation in mutations:
        if not mutation.residue.generic_number: continue  # cant map those without display numbers
        if mutation.residue.generic_number.label not in mutations_list: mutations_list[
            mutation.residue.generic_number.label] = []
        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = ''
        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = ''
        mutations_list[mutation.residue.generic_number.label].append(
            [mutation.foldchange, ligand.replace("'", "\\'"), qual])

    interaction_list = {} ###FIXME - always empty
    try:
        pc = ProteinConformation.objects.get(protein__family__slug=slug, protein__sequence_type__slug='consensus')
    except ProteinConformation.DoesNotExist:
        try:
            pc = ProteinConformation.objects.get(protein__family__slug=slug, protein__species_id=1,
                                                 protein__sequence_type__slug='wt')
        except:
            pc = None
            p = None
    p = pc.protein
    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    jsondata = {}
    jsondata_interaction = {}
    for r in residues:
        if r.generic_number:
            if r.generic_number.label in mutations_list:
                jsondata[r.sequence_number] = [mutations_list[r.generic_number.label]]
            if r.generic_number.label in interaction_list:
                jsondata_interaction[r.sequence_number] = interaction_list[r.generic_number.label]

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True

        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment)  # skip cells following title (which has colspan > 1)

        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []

        r_buffer.append((r, segment_title, title_cell_skip))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)
    
    context = {'pf': pf, 'families': families, 'structures': structures, 'no_of_proteins': no_of_proteins,
               'no_of_human_proteins': no_of_human_proteins, 'mutations': mutations, 'r_chunks': r_chunks,
               'chunk_size': chunk_size, 'p': p, 'complex_structures': complex_structures}

    return render(request,
                  'signprot/family_details.html',
                  context
    )

@cache_page(60 * 60 * 24 * 2)
def Ginterface(request, protein=None):
    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=protein).prefetch_related(
        'protein_segment', 'display_generic_number', 'generic_number')
    SnakePlot = DrawSnakePlot(
        residuelist, "Class A (Rhodopsin)", protein, nobuttons=1)

    # TEST
    gprotein_residues = Residue.objects.filter(protein_conformation__protein__entry_name='gnaz_human').prefetch_related(
        'protein_segment', 'display_generic_number', 'generic_number')
    gproteinplot = DrawGproteinPlot(
        gprotein_residues, "Gprotein", protein)

    crystal = Structure.objects.get(pdb_code__index="3SN6")
    aa_names = definitions.AMINO_ACID_GROUP_NAMES_OLD
    names_aa = dict(zip(aa_names.values(), aa_names.keys()))
    names_aa['Polar (S/T)'] = 'pol_short'
    names_aa['Polar (N/Q/H)'] = 'pol_long'

    residues_browser = [
        {'pos': 135, 'aa': 'I', 'gprotseg': "H5", 'segment': 'TM3', 'ligand': 'Gs', 'type': aa_names['hp'],
         'gpcrdb': '3.54x54', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},
        {'pos': 136, 'aa': 'T', 'gprotseg': "H5", 'segment': 'TM3', 'ligand': 'Gs', 'type': 'Polar (S/T)',
         'gpcrdb': '3.55x55', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},
        {'pos': 139, 'aa': 'F', 'gprotseg': "H5", 'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic',
         'gpcrdb': '34.51x51', 'gpnum': 'G.H5.8', 'gpaa': 'F376', 'availability': 'interacting'},
        {'pos': 139, 'aa': 'F', 'gprotseg': "S1", 'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic',
         'gpcrdb': '34.51x51', 'gpnum': 'G.S1.2', 'gpaa': 'H41', 'availability': 'interacting'},
        {'pos': 141, 'aa': 'Y', 'gprotseg': "H5", 'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic',
         'gpcrdb': '34.53x53', 'gpnum': 'G.H5.19', 'gpaa': 'H387', 'availability': 'interacting'},
        {'pos': 225, 'aa': 'E', 'gprotseg': "H5", 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge',
         'gpcrdb': '5.64x64', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},
        {'pos': 225, 'aa': 'E', 'gprotseg': "H5", 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge',
         'gpcrdb': '5.64x64', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},
        {'pos': 229, 'aa': 'Q', 'gprotseg': "H5", 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)',
         'gpcrdb': '5.68x68', 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'},
        {'pos': 229, 'aa': 'Q', 'gprotseg': "H5", 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)',
         'gpcrdb': '5.68x68', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},
        {'pos': 229, 'aa': 'Q', 'gprotseg': "H5", 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)',
         'gpcrdb': '5.68x68', 'gpnum': 'G.H5.17', 'gpaa': 'R385', 'availability': 'interacting'},
        {'pos': 274, 'aa': 'T', 'gprotseg': "H5", 'segment': 'TM6', 'ligand': 'Gs', 'type': 'Polar (S/T)',
         'gpcrdb': '6.36x36', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'},
        {'pos': 328, 'aa': 'R', 'gprotseg': "H5", 'segment': 'TM7', 'ligand': 'Gs', 'type': 'Positive charge',
         'gpcrdb': '7.55x55', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'},
        {'pos': 232, 'aa': 'K', 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Positive charge', 'gpcrdb': '5.71x71',
         'gprotseg': "H5", 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'}]

    # accessible_gn = ['3.50x50', '3.53x53', '3.54x54', '3.55x55', '34.50x50', '34.51x51', '34.53x53', '34.54x54', '5.61x61', '5.64x64', '5.65x65', '5.67x67', '5.68x68', '5.71x71', '5.72x72', '5.74x74', '5.75x75', '6.29x29', '6.32x32', '6.33x33', '6.36x36', '6.37x37', '7.55x55', '8.48x48', '8.49x49']

    accessible_gn = ['3.50x50', '3.53x53', '3.54x54', '3.55x55', '3.56x56', '34.50x50', '34.51x51', '34.52x52',
                     '34.53x53', '34.54x54', '34.55x55', '34.56x56', '34.57x57', '5.61x61', '5.64x64', '5.65x65',
                     '5.66x66', '5.67x67', '5.68x68', '5.69x69', '5.71x71', '5.72x72', '5.74x74', '5.75x75', '6.25x25',
                     '6.26x26', '6.28x28', '6.29x29', '6.32x32', '6.33x33', '6.36x36', '6.37x37', '6.40x40', '7.55x55',
                     '7.56x56', '8.47x47', '8.48x48', '8.49x49', '8.51x51']

    exchange_table = OrderedDict([('hp', ('V', 'I', 'L', 'M')),
                                  ('ar', ('F', 'H', 'W', 'Y')),
                                  ('pol_short', ('S', 'T')),  # Short/hydroxy
                                  ('pol_long', ('N', 'Q', 'H')),  # Amino-like (both donor and acceptor
                                  ('neg', ('D', 'E')),
                                  ('pos', ('K', 'R'))])

    interacting_gn = []

    accessible_pos = list(
        residuelist.filter(display_generic_number__label__in=accessible_gn).values_list('sequence_number', flat=True))

    # Which of the Gs interacting_pos are conserved?
    GS_none_equivalent_interacting_pos = []
    GS_none_equivalent_interacting_gn = []

    for interaction in residues_browser:
        interacting_gn.append(interaction['gpcrdb'])
        gs_b2_interaction_type_long = (
            next((item['type'] for item in residues_browser if item['gpcrdb'] == interaction['gpcrdb']), None))

        interacting_aa = residuelist.filter(display_generic_number__label__in=[interaction['gpcrdb']]).values_list(
            'amino_acid', flat=True)

        if interacting_aa:
            interaction['aa'] = interacting_aa[0]
            pos = \
                residuelist.filter(display_generic_number__label__in=[interaction['gpcrdb']]).values_list(
                    'sequence_number',
                    flat=True)[0]
            interaction['pos'] = pos

            feature = names_aa[gs_b2_interaction_type_long]

            if interacting_aa[0] not in exchange_table[feature]:
                GS_none_equivalent_interacting_pos.append(pos)
                GS_none_equivalent_interacting_gn.append(interaction['gpcrdb'])

    GS_equivalent_interacting_pos = list(
        residuelist.filter(display_generic_number__label__in=interacting_gn).values_list('sequence_number', flat=True))

    gProteinData = ProteinGProteinPair.objects.filter(protein__entry_name=protein)

    primary = []
    secondary = []

    for entry in gProteinData:
        if entry.transduction == 'primary':
            primary.append((entry.g_protein.name.replace("Gs", "G<sub>s</sub>").replace("Gi", "G<sub>i</sub>").replace(
                "Go", "G<sub>o</sub>").replace("G11", "G<sub>11</sub>").replace("G12", "G<sub>12</sub>").replace("G13",
                                                                                                                 "G<sub>13</sub>").replace(
                "Gq", "G<sub>q</sub>").replace("G", "G&alpha;"), entry.g_protein.slug))
        elif entry.transduction == 'secondary':
            secondary.append((
                entry.g_protein.name.replace("Gs", "G<sub>s</sub>").replace("Gi", "G<sub>i</sub>").replace(
                    "Go", "G<sub>o</sub>").replace("G11", "G<sub>11</sub>").replace("G12",
                                                                                    "G<sub>12</sub>").replace(
                    "G13", "G<sub>13</sub>").replace("Gq", "G<sub>q</sub>").replace("G", "G&alpha;"),
                entry.g_protein.slug))

    return render(request,
                  'signprot/ginterface.html',
                  {'pdbname': '3SN6',
                   'snakeplot': SnakePlot,
                   'gproteinplot': gproteinplot,
                   'crystal': crystal,
                   'interacting_equivalent': GS_equivalent_interacting_pos,
                   'interacting_none_equivalent': GS_none_equivalent_interacting_pos,
                   'accessible': accessible_pos,
                   'residues': residues_browser,
                   'mapped_protein': protein,
                   'interacting_gn': GS_none_equivalent_interacting_gn,
                   'primary_Gprotein': set(primary),
                   'secondary_Gprotein': set(secondary)}
    )

def ajaxInterface(request, slug, **response_kwargs):
    name_of_cache = 'ajaxInterface_' + slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:

        p = Protein.objects.filter(entry_name=slug).get()

        if p.family.slug.startswith('200'):
            rsets = ResiduePositionSet.objects.get(name="Arrestin interface")
        else:
            rsets = ResiduePositionSet.objects.get(name="Gprotein Barcode")

        jsondata = {}
        for x, residue in enumerate(rsets.residue_position.all()):
            try:
                pos = str(list(Residue.objects.filter(protein_conformation__protein__entry_name=slug,
                                                      display_generic_number__label=residue.label))[0])
            except:
                print("Protein has no residue position at", residue.label)
            a = pos[1:]

            jsondata[a] = [5, 'Receptor interface position', residue.label]

        jsondata = json.dumps(jsondata)

    cache.set(name_of_cache, jsondata, 60 * 60 * 24 * 2)  # two days timeout on cache

    response_kwargs['content_type'] = 'application/json'

    return HttpResponse(jsondata, **response_kwargs)

def ajaxBarcode(request, slug, cutoff, **response_kwargs):
    name_of_cache = 'ajaxBarcode_' + slug + cutoff

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        selectivity_pos = list(
            SignprotBarcode.objects.filter(protein__entry_name=slug, seq_identity__gte=cutoff).values_list(
                'residue__display_generic_number__label', flat=True))

        conserved = list(SignprotBarcode.objects.filter(protein__entry_name=slug, paralog_score__gte=cutoff,
                                                        seq_identity__gte=cutoff).prefetch_related(
            'residue__display_generic_number').values_list('residue__display_generic_number__label', flat=True))

        na_data = list(
            SignprotBarcode.objects.filter(protein__entry_name=slug, seq_identity=0, paralog_score=0).values_list(
                'residue__display_generic_number__label', flat=True))

        all_positions = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related(
            'display_generic_number')

        for res in all_positions:
            cgn = str(res.generic_number)
            res = str(res.sequence_number)
            if cgn in conserved:
                jsondata[res] = [0, 'Conserved', cgn]
            elif cgn in selectivity_pos and cgn not in conserved:
                jsondata[res] = [1, 'Selectivity determining', cgn]
            elif cgn in na_data:
                jsondata[res] = [3, 'NA', cgn]
            else:
                jsondata[res] = [2, 'Evolutionary neutral', cgn]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 60 * 60 * 24 * 2)  # two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24 * 2)
def StructureInfo(request, pdbname):
    """
    Show structure details
    """
    protein = Protein.objects.get(signprotstructure__pdb_code__index=pdbname)

    crystal = SignprotStructure.objects.get(pdb_code__index=pdbname)

    return render(request,
                  'signprot/structure_info.html',
                  {'pdbname': pdbname,
                   'protein': protein,
                   'crystal': crystal}
    )

# @cache_page(60*60*24*2)
def signprotdetail(request, slug):
    # get protein

    slug = slug.lower()
    p = Protein.objects.prefetch_related('web_links__web_resource').get(entry_name=slug, sequence_type__slug='wt')

    # get family list
    pf = p.family
    families = [pf.name]
    while pf.parent.parent:
        families.append(pf.parent.name)
        pf = pf.parent
    families.reverse()

    # get protein aliases
    aliases = ProteinAlias.objects.filter(protein=p).values_list('name', flat=True)

    # get genes
    genes = Gene.objects.filter(proteins=p).values_list('name', flat=True)
    gene = genes[0]
    alt_genes = genes[1:]

    # get structures of this signal protein
    structures = SignprotStructure.objects.filter(protein=p)
    complex_structures = SignprotComplex.objects.filter(protein=p)
    
    # mutations
    mutations = MutationExperiment.objects.filter(protein=p)

    # get residues
    pc = ProteinConformation.objects.get(protein=p)

    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True

        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment)  # skip cells following title (which has colspan > 1)

        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []

        r_buffer.append((r, segment_title, title_cell_skip))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)
    context = {'p': p, 'families': families, 'r_chunks': r_chunks, 'chunk_size': chunk_size, 'aliases': aliases,
               'gene': gene, 'alt_genes': alt_genes, 'structures': structures, 'complex_structures': complex_structures,
               'mutations': mutations}
    
    return render(request,
                  'signprot/signprot_details.html',
                  context
    )

def sort_a_by_b(a, b, remove_invalid=False):
    '''Sort one list based on the order of elements from another list'''
    # https://stackoverflow.com/q/12814667
    # a = ['alpha_mock', 'van-der-waals', 'ionic']
    # b = ['ionic', 'aromatic', 'hydrophobic', 'polar', 'van-der-waals', 'alpha_mock']
    # sort_a_by_b(a,b) -> ['ionic', 'van-der-waals', 'alpha_mock']
    if remove_invalid:
        a = [a_elem for a_elem in a if a_elem in b]
    return sorted(a, key=lambda x: b.index(x))

def interface_dataset():
    # correct receptor entry names - the ones with '_a' appended
    complex_objs = SignprotComplex.objects.prefetch_related('structure__protein_conformation__protein')

    # TOFIX: Current workaround is forcing _a to pdb for indicating alpha-subunit
    #complex_names = [complex_obj.structure.protein_conformation.protein.entry_name + '_' + complex_obj.alpha.lower() for
    #                 complex_obj in complex_objs]
    complex_names = [complex_obj.structure.protein_conformation.protein.entry_name + '_a' for
                     complex_obj in complex_objs]

    complex_struc_ids = [co.structure_id for co in complex_objs]
    # protein conformations for those
    prot_conf = ProteinConformation.objects.filter(protein__entry_name__in=complex_names).values_list('id', flat=True)

    interaction_sort_order = [
        "ionic",
        "aromatic",
        "polar",
        "hydrophobic",
        "van-der-waals",
    ]

    # getting all the signal protein residues for those protein conformations
    prot_residues = Residue.objects.filter(
        protein_conformation__in=prot_conf
    ).values_list('id', flat=True)

    interactions = InteractingResiduePair.objects.filter(
        Q(res1__in=prot_residues) | Q(res2__in=prot_residues),
        referenced_structure__in=complex_struc_ids
    ).exclude(
        Q(res1__in=prot_residues) & Q(res2__in=prot_residues)
    ).prefetch_related(
        'interaction__interaction_type',
        'referenced_structure__pdb_code__index',
        'referenced_structure__signprot_complex__protein__entry_name',
        'referenced_structure__protein_conformation__protein__parent__entry_name',
        'res1__amino_acid',
        'res1__sequence_number',
        'res1__generic_number__label',
        'res2__amino_acid',
        'res2__sequence_number',
        'res2__generic_number__label',
    ).order_by(
        'res1__generic_number__label',
        'res2__generic_number__label'
    ).values(
        int_id=F('id'),
        int_ty=ArrayAgg(
            'interaction__interaction_type',
            distinct=True,
            # ordering=interaction_sort_order
        ),

        pdb_id=F('referenced_structure__pdb_code__index'),
        conf_id=F('referenced_structure__protein_conformation_id'),
        gprot=F('referenced_structure__signprot_complex__protein__entry_name'),
        entry_name=F('referenced_structure__protein_conformation__protein__parent__entry_name'),

        rec_aa=F('res1__amino_acid'),
        rec_pos=F('res1__sequence_number'),
        rec_gn=F('res1__generic_number__label'),

        sig_aa=F('res2__amino_acid'),
        sig_pos=F('res2__sequence_number'),
        sig_gn=F('res2__generic_number__label')
    )

    conf_ids = set()
    for i in interactions:
        i['int_ty'] = sort_a_by_b(i['int_ty'], interaction_sort_order)
        conf_ids.update([i['conf_id']])

    return list(conf_ids), list(interactions)

# @cache_page(60*60*24*2)
def InteractionMatrix(request):
    prot_conf_ids, dataset = interface_dataset()

    gprotein_order = ProteinSegment.objects.filter(proteinfamily='Alpha').values('id', 'slug')
    receptor_order = ['N', '1', '12', '2', '23', '3', '34', '4', '45', '5', '56', '6', '67', '7', '78', '8', 'C']

    struc = SignprotComplex.objects.prefetch_related(
        'structure__pdb_code',
        'structure__stabilizing_agents',
        'structure__protein_conformation__protein__species',
        'structure__protein_conformation__protein__parent__parent__parent',
        'structure__protein_conformation__protein__family__parent__parent__parent__parent',
        'structure__stabilizing_agents',
        'structure__signprot_complex__protein__family__parent__parent__parent__parent',
    )

    complex_info = []
    for s in struc:
        r = {}
        s = s.structure
        r['pdb_id'] = s.pdb_code.index
        r['name'] = s.protein_conformation.protein.parent.short()
        r['entry_name'] = s.protein_conformation.protein.parent.entry_name
        r['class'] = s.protein_conformation.protein.get_protein_class()
        r['family'] = s.protein_conformation.protein.get_protein_family()
        r['conf_id'] = s.protein_conformation.id
        r['organism'] = s.protein_conformation.protein.species.common_name
        try:
            r['gprot'] = s.get_stab_agents_gproteins()
        except Exception:
            r['gprot'] = ''
        try:
            r['gprot_class'] = s.get_signprot_gprot_family()
        except Exception:
            r['gprot_class'] = ''
        complex_info.append(r)

    remaining_residues = Residue.objects.filter(
        protein_conformation_id__in=prot_conf_ids,
    ).prefetch_related(
        "protein_conformation",
        "protein_conformation__protein",
        "protein_conformation__structure"
    ).values(
        rec_id=F('protein_conformation__protein__id'),
        name=F('protein_conformation__protein__parent__name'),
        entry_name=F('protein_conformation__protein__parent__entry_name'),
        pdb_id=F('protein_conformation__structure__pdb_code__index'),
        rec_aa=F('amino_acid'),
        rec_gn=F('generic_number__label'),
    ).exclude(
        Q(rec_gn=None)
    )

    context = {
        'interactions': json.dumps(dataset),
        'interactions_metadata': json.dumps(complex_info),
        'non_interactions': json.dumps(list(remaining_residues)),
        'gprot': json.dumps(list(gprotein_order)),
        'receptor': json.dumps(receptor_order),
    }

    request.session['signature'] = None
    request.session.modified = True

    return render(request,
                  'signprot/matrix.html',
                  context
    )

def IMSequenceSignature(request):
    """Accept set of proteins + generic numbers and calculate the signature for those"""
    t1 = time.time()

    pos_set_in = get_entry_names(request)
    ignore_in_alignment = get_ignore_info(request)
    segments = get_protein_segments(request)
    if len(segments) == 0:
        segments = list(ResidueGenericNumberEquivalent.objects.filter(scheme__slug__in=['gpcrdba']))

    # get pos objects
    pos_set = Protein.objects.filter(entry_name__in=pos_set_in).select_related('residue_numbering_scheme', 'species')

    # Calculate Sequence Signature
    signature = SequenceSignature()

    signature.setup_alignments_signprot(segments, pos_set, ignore_in_alignment=ignore_in_alignment)
    signature.calculate_signature_onesided()
    # preprocess data for return
    signature_data = signature.prepare_display_data_onesided()

    # FEATURES AND REGIONS
    feats = [feature for feature in signature_data['a_pos'].features_combo]

    # GET GENERIC NUMBERS
    generic_numbers = get_generic_numbers(signature_data)

    # FEATURE FREQUENCIES
    signature_features = get_signature_features(signature_data, generic_numbers, feats)
    grouped_features = group_signature_features(signature_features)

    # # FEATURE CONSENSUS
    # generic_numbers_flat = list(chain.from_iterable(generic_numbers))
    # sigcons = get_signature_consensus(signature_data, generic_numbers_flat)

    # rec_class = pos_set[0].get_protein_class()

    # dump = {
    #     'rec_class': rec_class,
    #     'signature': signature,
    #     'consensus': signature_data,
    #     }
    # with open('signprot/notebooks/interface_pickles/{}.p'.format(rec_class), 'wb+') as out_file:
    #     pickle.dump(dump, out_file)

    # pass back to front
    res = {
        # 'cons': sigcons,
        'feat_ungrouped': signature_features,
        'feat': grouped_features,
    }

    request.session['signature'] = signature.prepare_session_data()
    request.session.modified = True

    t2 = time.time()
    print('Runtime: {}'.format((t2 - t1) * 1000.0))

    return JsonResponse(res, safe=False)

def IMSignatureMatch(request):
    '''Take the signature stored in the session and query the db'''
    signature_data = request.session.get('signature')
    ss_pos = get_entry_names(request)
    cutoff = request.POST.get('cutoff')

    request.session['ss_pos'] = ss_pos
    request.session['cutoff'] = cutoff

    pos_set = Protein.objects.filter(entry_name__in=ss_pos).select_related('residue_numbering_scheme', 'species')
    pos_set = [protein for protein in pos_set]
    pfam = [protein.family.slug[:3] for protein in pos_set]

    signature_match = SignatureMatch(
        signature_data['common_positions'],
        signature_data['numbering_schemes'],
        signature_data['common_segments'],
        signature_data['diff_matrix'],
        pos_set,
        # pos_set,
        cutoff=0,
        signprot=True
    )

    maj_pfam = Counter(pfam).most_common()[0][0]
    signature_match.score_protein_class(maj_pfam, signprot=True)
    # request.session['signature_match'] = signature_match

    signature_match = {
        'scores': signature_match.protein_report,
        'scores_pos': signature_match.scores_pos,
        # 'scores_neg': signature_match.scores_neg,
        'protein_signatures': signature_match.protein_signatures,
        'signatures_pos': signature_match.signatures_pos,
        # 'signatures_neg': signature_match.signatures_neg,
        'signature_filtered': signature_match.signature_consensus,
        'relevant_gn': signature_match.relevant_gn,
        'relevant_segments': signature_match.relevant_segments,
        'numbering_schemes': signature_match.schemes,
    }

    signature_match = prepare_signature_match(signature_match)
    return JsonResponse(signature_match, safe=False)

def render_IMSigMat(request):
    # signature_match = request.session.get('signature_match')
    signature_data = request.session.get('signature')
    ss_pos = request.session.get('ss_pos')
    cutoff = request.session.get('cutoff')

    pos_set = Protein.objects.filter(entry_name__in=ss_pos).select_related('residue_numbering_scheme', 'species')
    pos_set = [protein for protein in pos_set]
    pfam = [protein.family.slug[:3] for protein in pos_set]

    signature_match = SignatureMatch(
        signature_data['common_positions'],
        signature_data['numbering_schemes'],
        signature_data['common_segments'],
        signature_data['diff_matrix'],
        pos_set,
        # pos_set,
        cutoff=0,
        signprot=True
    )

    maj_pfam = Counter(pfam).most_common()[0][0]
    signature_match.score_protein_class(maj_pfam, signprot=True)

    response = render(
        request,
        'signprot/signature_match.html',
        {'scores': signature_match}
    )
    return response
