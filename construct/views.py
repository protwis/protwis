from django.shortcuts import render
from django.views.generic import TemplateView, View
from django.http import HttpResponse
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings


from common.diagrams_gpcr import DrawSnakePlot
from common.views import AbsTargetSelection
from common.definitions import AMINO_ACIDS, AMINO_ACID_GROUPS
from construct.models import *
from construct.functions import *
from protein.models import Protein, ProteinConformation,ProteinSegment
from structure.models import Structure
from mutation.models import Mutation
from residue.models import ResiduePositionSet
from construct.tool import *


from datetime import datetime
import json
import copy
from collections import OrderedDict


Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')


# Create your views here.
#@cache_page(60 * 60 * 24)
def detail(request, slug):

    # get constructs
    c = Construct.objects.filter(name=slug).prefetch_related(
        'modifications', 'mutations', 'deletions','insertions','insertions__insert_type',
        'expression','solubilization','purification','crystallization','crystal').all()[0]

    # get residues
    residues = Residue.objects.filter(protein_conformation__protein=c.protein).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    residues_lookup = {}
    for r in residues:
        residues_lookup[r.sequence_number] = r

    schematics = c.schematic()

    snake = c.snake()

    # snake = cache.get(c.name+'_snake')
    # snake= None
    # if snake==None:
    #     print(c.name+'_snake no cache')
    #     snake = cache.set(c.name+'_snake', DrawSnakePlot(residues,c.protein.get_protein_class(),str(c.protein),nobuttons = True), 60*60*24*2) #two days
    #     snake = cache.get(c.name+'_snake')
    # else:
    #     print(c.name+'_snake used cache')

    chunk_size = 10
    context = {'c':c, 'chunk_size': chunk_size, 'snake': snake, 'annotations': json.dumps(schematics['annotations']), 'schematics': schematics, 'residues_lookup': residues_lookup}
    return render(request,'construct/construct_detail.html',context)

class ConstructStatistics(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct/statistics.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructStatistics, self).get_context_data(**kwargs)
        cons = Construct.objects.all().prefetch_related(
            "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
            "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

        #PREPARE DATA
        proteins = Construct.objects.all().values_list('protein', flat = True)
        pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['1x50','8x50','5x50','6x50']).values_list('protein__entry_name','residue__sequence_number','residue__generic_number__label')
        #print(pconfs)
        x50s = {}
        for pc in pconfs:
            if pc[0] not in x50s:
                x50s[pc[0]] = {}
            x50s[pc[0]][pc[2]] = pc[1]

        #GRAB RESIDUES for mutations
        mutations = []
        positions = []
        proteins = []
        for c in cons:
            p = c.protein
            entry_name = p.entry_name
            p_class = p.family.slug.split('_')[0]
            pdb = c.crystal.pdb_code
            for mutation in c.mutations.all():
                if p.entry_name not in proteins:
                    proteins.append(entry_name)
                mutations.append((mutation,entry_name,pdb,p_class))
                if mutation.sequence_number not in positions:
                    positions.append(mutation.sequence_number)
        rs = Residue.objects.filter(protein_conformation__protein__entry_name__in=proteins, sequence_number__in=positions).prefetch_related('generic_number','protein_conformation__protein')

        rs_lookup = {}
        gns = []
        for r in rs:
            if not r.generic_number:
                continue #skip non GN
            entry_name = r.protein_conformation.protein.entry_name
            pos = r.sequence_number
            if entry_name not in rs_lookup:
                rs_lookup[entry_name] = {}
            if pos not in rs_lookup[entry_name]:
                rs_lookup[entry_name][pos] = r

        truncations = {}
        truncations['nterm'] = {}
        truncations['nterm_fusion'] = {}
        truncations['cterm'] = {}
        truncations['icl3'] = {}
        truncations['icl3_fusion'] = {}
        class_names = {}
        for c in cons:
            p = c.protein
            entry_name = p.entry_name
            p_class = p.family.slug.split('_')[0]
            if p_class not in class_names:
                class_names[p_class] = p.family.parent.parent.parent.name
            p_class_name = class_names[p_class]
            fusion_n = False
            fusion_icl3 = False

            fusion_position, fusions = c.fusion()

            for deletion in c.deletions.all():
                if deletion.end <= x50s[entry_name]['1x50']:
                    bw = "1."+str(50-x50s[entry_name]['1x50']+deletion.end)
                    #bw = bw + " " + str(x50s[entry_name]['1x50']-deletion.end)

                    position = 'nterm'
                    if fusion_position=='nterm':
                        position = 'nterm_fusion'

                    if p_class_name not in truncations[position]:
                        truncations[position][p_class_name] = {}
                    if bw not in truncations[position][p_class_name]:
                        truncations[position][p_class_name][bw] = []
                    if entry_name not in truncations[position][p_class_name][bw]:
                        truncations[position][p_class_name][bw].append(entry_name)
                if deletion.start >= x50s[entry_name]['8x50']:
                    bw = x50s[entry_name]['8x50']-deletion.start
                    bw = "8."+str(50-x50s[entry_name]['8x50']+deletion.start)
                    if p_class_name not in truncations['cterm']:
                        truncations['cterm'][p_class_name] = {}
                    if bw not in truncations['cterm'][p_class_name]:
                        truncations['cterm'][p_class_name][bw] = []
                    if entry_name not in truncations['cterm'][p_class_name][bw]:
                        truncations['cterm'][p_class_name][bw].append(entry_name)
                if deletion.start > x50s[entry_name]['5x50'] and deletion.start < x50s[entry_name]['6x50']:
                    bw = x50s[entry_name]['5x50']-deletion.start
                    bw = "5."+str(50-x50s[entry_name]['5x50']+deletion.start)
                    bw2 = "6."+str(50-x50s[entry_name]['6x50']+deletion.end)
                    bw_combine = bw+"-"+bw2
                    position = 'icl3'
                    if fusion_position=='icl3':
                        position = 'icl3_fusion'
                    if p_class_name not in truncations[position]:
                        truncations[position][p_class_name] = {}
                    if bw_combine not in truncations[position][p_class_name]:
                        truncations[position][p_class_name][bw_combine] = []
                    if entry_name not in truncations[position][p_class_name][bw_combine]:
                        truncations[position][p_class_name][bw_combine].append(entry_name)
        #print(truncations)
        #truncations = OrderedDict(truncations)
        ordered_truncations = OrderedDict()
        for segment, s_vals in sorted(truncations.items()):
            #print(segment)
            ordered_truncations[segment] = OrderedDict()
            for p_class, c_vals in sorted(s_vals.items()):
                #print(p_class) 
                ordered_truncations[segment][p_class] = OrderedDict()
                for pos, p_vals in sorted(c_vals.items(),key=lambda x: (len(x[1]),x[0]), reverse=True):
                    #print(pos, len(p_vals))
                    ordered_truncations[segment][p_class][pos] = p_vals

        #truncations =  OrderedDict(sorted(truncations.items(), key=lambda x: x[1]['hits'],reverse=True))
        #print(ordered_truncations)
        context['truncations'] = ordered_truncations

        mutation_list = OrderedDict()
        mutation_type = OrderedDict()
        mutation_wt = OrderedDict()
        mutation_mut = OrderedDict()
        for mutation in mutations:
            wt = mutation[0].wild_type_amino_acid
            mut = mutation[0].mutated_amino_acid
            entry_name = mutation[1]
            pos = mutation[0].sequence_number
            p_class = mutation[3]
            p_class = class_names[p_class]
            pdb = mutation[2]

            print(entry_name,"\t", pdb,"\t", pos,"\t", wt,"\t", mut)


            if p_class not in mutation_type:
                mutation_type[p_class] = OrderedDict()
            if wt+"=>"+mut not in mutation_type[p_class]:
                mutation_type[p_class][wt+"=>"+mut] = {'hits':0, 'proteins':[]}
            if entry_name not in mutation_type[p_class][wt+"=>"+mut]['proteins']:
                mutation_type[p_class][wt+"=>"+mut]['proteins'].append(entry_name)
                mutation_type[p_class][wt+"=>"+mut]['hits'] += 1


            if p_class not in mutation_wt:
                mutation_wt[p_class] = OrderedDict()
            if wt not in mutation_wt[p_class]:
                mutation_wt[p_class][wt] = {'hits':0, 'proteins':[]}
            if entry_name not in mutation_wt[p_class][wt]['proteins']:
                mutation_wt[p_class][wt]['proteins'].append(entry_name)
                mutation_wt[p_class][wt]['hits'] += 1

            if p_class not in mutation_mut:
                mutation_mut[p_class] = OrderedDict()
            if mut not in mutation_mut[p_class]:
                mutation_mut[p_class][mut] = {'hits':0, 'proteins':[]}
            if entry_name not in mutation_mut[p_class][mut]['proteins']:
                mutation_mut[p_class][mut]['proteins'].append(entry_name)
                mutation_mut[p_class][mut]['hits'] += 1

            if entry_name not in rs_lookup:
                continue
            if pos not in rs_lookup[entry_name]:
                continue
            gn = rs_lookup[entry_name][pos].generic_number.label

            if p_class not in mutation_list:
                mutation_list[p_class] = OrderedDict()
            if gn not in mutation_list[p_class]:
                mutation_list[p_class][gn] = {'proteins':[], 'hits':0, 'mutation':[]}
            if entry_name not in mutation_list[p_class][gn]['proteins']:
                mutation_list[p_class][gn]['proteins'].append(entry_name)
                mutation_list[p_class][gn]['hits'] += 1  
                mutation_list[p_class][gn]['mutation'].append((mutation[0].wild_type_amino_acid,mutation[0].mutated_amino_acid))



        for p_class, values in mutation_list.items():
            for gn, vals in values.items():
                if vals['hits']<2:
                    pass
                    #values.pop(gn, None)
            mutation_list[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))
            #mutation_list = OrderedDict(sorted(mutation_list.items(), key=lambda x: x[1]['hits'],reverse=True))

        for p_class, values in mutation_type.items():
            mutation_type[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))

        for p_class, values in mutation_wt.items():
            mutation_wt[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))

        for p_class, values in mutation_mut.items():
            mutation_mut[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))

        context['mutation_list'] = mutation_list
        context['mutation_type'] = mutation_type
        context['mutation_wt'] = mutation_wt
        context['mutation_mut'] = mutation_mut

        for c in cons:
            pass

        return context


class ConstructTable(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct/residuetable.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructTable, self).get_context_data(**kwargs)
        cons = Construct.objects.all().prefetch_related(
            "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
            "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

        #PREPARE DATA
        proteins = Construct.objects.all().values_list('protein', flat = True)

        #GRAB RESIDUES for mutations
        mutations = []
        positions = []
        proteins = []
        class_names = {}
        classes = []
        for c in cons:
            p = c.protein
            entry_name = p.entry_name
            p_class = p.family.slug.split('_')[0]
            if p_class not in classes:
                classes.append(p_class)
            pdb = c.crystal.pdb_code
            for mutation in c.mutations.all():
                if p.entry_name not in proteins:
                    proteins.append(entry_name)
                mutations.append((mutation,entry_name,pdb,p_class,c.name))
                if mutation.sequence_number not in positions:
                    positions.append(mutation.sequence_number)
        rs = Residue.objects.filter(protein_conformation__protein__entry_name__in=proteins, sequence_number__in=positions).prefetch_related('generic_number','protein_conformation__protein','protein_segment')

        #print(classes)

        excluded_segment = ['C-term','N-term']
        list(settings.REFERENCE_POSITIONS.keys())
        align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()

        amino_acids_stats = {}
        amino_acids_groups_stats = {}
        accessible_in_class = {}
        for c in classes:

            #if c=='001':
            #    continue
            amino_acids_stats[c] = {}
            amino_acids_groups_stats[c] = {}
            accessible_in_class[c] = []

            if c =='001':
                residue_set_name = 'Class A binding pocket'
            elif c=='004':
                residue_set_name = 'Class C binding pocket'
            elif c=='005':
                residue_set_name = 'Class F binding pocket'
            else:
                residue_set_name = ''

            if residue_set_name:
                rset = ResiduePositionSet.objects.get(name=residue_set_name)
                for residue in rset.residue_position.all():
                    accessible_in_class[c].append(residue.label)
            #print(accessible_in_class)

            alignment_proteins = Protein.objects.filter(family__slug__startswith=c, species__common_name='Human', source__name='SWISSPROT')
            #print(c,len(alignment_proteins))

            a = Alignment()

            a.load_proteins(alignment_proteins)

            a.load_segments(align_segments) #get all segments to make correct diagrams

            # build the alignment data matrix
            a.build_alignment()

            # calculate consensus sequence + amino acid and feature frequency
            a.calculate_statistics()

            #print(a.amino_acid_stats)
                        # {% for ns, segments in a.generic_numbers.items %}
                        # <tr>
                        #     {% for s, num in segments.items %}
                        #         {% for n, dn in num.items %}
                        #             {% if 'Common G-alpha numbering scheme' in a.numbering_schemes.0 %}
                        #             <td class="ali-td-generic-num">{{ dn|make_list|slice:'2:'|join:''}}</td>
                        #             {% else %}
                        #             <td class="ali-td-generic-num">{{ dn|safe }}</td>
                        #             {% endif %}
                        #         {% endfor %}
                        #         <td class="ali-td">&nbsp;</td>
                        #     {% endfor %}
                        # </tr>
                        # {% endfor %}
            s_id = 0
            a_id = 0
            for ns, segments in a.generic_numbers.items():
                for s, num in segments.items():
                    for n, dn in num.items():
                        temp = []
                        temp2 = []
                        for i, aa in enumerate(AMINO_ACIDS):
                            temp.append(a.amino_acid_stats[i][s_id][a_id])

                        for i, aa in enumerate(AMINO_ACID_GROUPS):
                            temp2.append(a.feature_stats[i][s_id][a_id])
                        amino_acids_stats[c][n] = temp
                        amino_acids_groups_stats[c][n] = temp2
                    a_id += 1
                s_id += 1

        rs_lookup = {}
        gns = []
        for r in rs:
            entry_name = r.protein_conformation.protein.entry_name
            pos = r.sequence_number
            segment = r.protein_segment.slug
            if entry_name not in rs_lookup:
                rs_lookup[entry_name] = {}
            if pos not in rs_lookup[entry_name]:
                rs_lookup[entry_name][pos] = r

        count_per_gn = {}
        for mutation in mutations:
            entry_name = mutation[1]
            pos = mutation[0].sequence_number
            if rs_lookup[entry_name][pos].generic_number:
                gn = rs_lookup[entry_name][pos].generic_number.label
                if gn not in count_per_gn:
                    count_per_gn[gn] = {'hits': 0, 'proteins': []}
                if entry_name not in count_per_gn[gn]['proteins']:
                    count_per_gn[gn]['proteins'].append(entry_name)
                    count_per_gn[gn]['hits'] += 1

        #print(count_per_gn)

        mutation_list = []
        for mutation in mutations:
            wt = mutation[0].wild_type_amino_acid
            mut = mutation[0].mutated_amino_acid
            entry_name = mutation[1]
            pdb = mutation[2]
            cname = mutation[4]
            pos = mutation[0].sequence_number
            p_class = mutation[3]
            if p_class not in class_names:
                class_names[p_class] = p.family.parent.parent.parent.name
            p_class_name = class_names[p_class]
            p_class = class_names[p_class]



            if entry_name not in rs_lookup:
                continue
            if pos not in rs_lookup[entry_name]:
                continue
            segment = rs_lookup[entry_name][pos].protein_segment.slug
            if rs_lookup[entry_name][pos].generic_number:
                gn = rs_lookup[entry_name][pos].generic_number.label
                stats = amino_acids_stats[mutation[3]][gn]
                stats2 = amino_acids_groups_stats[mutation[3]][gn]
                if gn in accessible_in_class[mutation[3]]:
                    accessible = 'yes'
                else:
                    accessible = 'no'
                count = count_per_gn[gn]['hits']
            else:
                gn = ''
                stats = ''
                stats2 = ''
                accessible = 'N/A'



            mutation_list.append({'entry_name':entry_name,'pdb':pdb,'cname':cname, 'segment':segment,'pos': pos, 'gn': gn, 'wt': wt, 'mut': mut,'p_class': p_class, 'amino_acids':stats, 'amino_acids_groups':stats2, 'accessible': accessible, 'count': count})


        context['amino_acids'] = AMINO_ACIDS
        context['amino_groups'] = AMINO_ACID_GROUPS
        context['mutation_list'] = mutation_list
        #print(mutation_list)

        return context

class ConstructMutations(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct/mutations.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructMutations, self).get_context_data(**kwargs)
        cons = Construct.objects.all().prefetch_related(
            "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
            "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

        #PREPARE DATA
        proteins = Construct.objects.all().values_list('protein', flat = True)

        #GRAB RESIDUES for mutations
        mutations = []
        positions = []
        proteins = []
        class_names = {}
        for c in cons:
            p = c.protein
            entry_name = p.entry_name
            p_class = p.family.slug.split('_')[0]
            pdb = c.crystal.pdb_code
            for mutation in c.mutations.all():
                if p.entry_name not in proteins:
                    proteins.append(entry_name)
                mutations.append((mutation,entry_name,pdb,p_class,c.name))
                if mutation.sequence_number not in positions:
                    positions.append(mutation.sequence_number)
        rs = Residue.objects.filter(protein_conformation__protein__entry_name__in=proteins, sequence_number__in=positions).prefetch_related('generic_number','protein_conformation__protein','protein_segment')

        rs_lookup = {}
        gns = []
        for r in rs:
            entry_name = r.protein_conformation.protein.entry_name
            pos = r.sequence_number
            segment = r.protein_segment.slug
            if entry_name not in rs_lookup:
                rs_lookup[entry_name] = {}
            if pos not in rs_lookup[entry_name]:
                rs_lookup[entry_name][pos] = r

        mutation_list = []
        for mutation in mutations:
            wt = mutation[0].wild_type_amino_acid
            mut = mutation[0].mutated_amino_acid
            entry_name = mutation[1]
            pdb = mutation[2]
            cname = mutation[4]
            pos = mutation[0].sequence_number
            p_class = mutation[3]
            if p_class not in class_names:
                class_names[p_class] = p.family.parent.parent.parent.name
            p_class_name = class_names[p_class]
            p_class = class_names[p_class]


            if entry_name not in rs_lookup:
                continue
            if pos not in rs_lookup[entry_name]:
                continue
            segment = rs_lookup[entry_name][pos].protein_segment.slug
            if rs_lookup[entry_name][pos].generic_number:
                gn = rs_lookup[entry_name][pos].generic_number.label
            else:
                gn = ''


            mutation_list.append({'entry_name':entry_name,'pdb':pdb,'cname':cname, 'segment':segment,'pos': pos, 'gn': gn, 'wt': wt, 'mut': mut,'p_class': p_class})


        context['mutation_list'] = mutation_list
        print(mutation_list)

        return context

def fetch_all_pdb(request):

    structures = Structure.objects.all()

    for s in structures:
        pdbname = str(s)
        print(pdbname)
        failed = []
        try:
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            d = fetch_pdb_info(pdbname,protein)

            #delete before adding new
            Construct.objects.filter(name=d['construct_crystal']['pdb_name']).delete()
            add_construct(d)
        except:
            print(pdbname,'failed')
            failed.append(pdbname)


    context = {'failed':failed}

    return render(request,'pdb_all.html',context)

def fetch_pdb(request, slug):

    protein = Protein.objects.filter(entry_name=slug.lower()).get()

    d = fetch_pdb_info(slug,protein)

    #delete before adding new
    print(d['construct_crystal']['pdb_name'])
    Construct.objects.filter(name__iexact=d['construct_crystal']['pdb_name']).delete()
    add_construct(d)

    cache.delete(d['construct_crystal']['pdb_name']+'_schematics')
    cache.delete(d['construct_crystal']['pdb_name']+'_snake')
    context = {'d':d}

    return render(request,'pdb_fetch.html',context)

def fetch_pdb_for_webform(request, slug, **response_kwargs):

    slug = slug.lower()
    protein = Protein.objects.filter(entry_name=slug).get()

    d = fetch_pdb_info(slug,protein)
    d = convert_ordered_to_disordered_annotation(d)

    jsondata = json.dumps(d)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

class ConstructBrowser(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct_browser.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructBrowser, self).get_context_data(**kwargs)
        try:
            cons = Construct.objects.all().prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

            #context['constructs'] = cache.get('construct_browser')
            #if context['constructs']==None:
            context['constructs'] = []
            for c in cons:
                c.schematics = c.schematic()
                context['constructs'].append(c)

            #cache.set('construct_browser', context['constructs'], 60*60*24*2) #two days
            # else:
            #     print('construct_browser used cache')

        except Construct.DoesNotExist as e:
            pass

        return context

class design(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 1
    # docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Get construct suggestions based on published constructs.'

    # Middle section
    numbering_schemes = False
    filters = False
    search = True
    title = "Select a receptor"

    template_name = 'designselection.html'
    type_of_selection = 'targets'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # Buttons
    buttons = {
        'continue': {
            'label': 'Show results',
            'onclick': 'submitupload()',
            'color': 'success',
            #'url': 'calculate/'
        }
    }

    redirect_on_select = False

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        return context

@csrf_exempt #jquery send post, so no csrf
def align(request):

    ids = json.loads(request.POST.get('ids'))

    c_ids = []
    s_ids = []
    for i in ids:
        if i.startswith('align'):
            s_ids.append(i.split('_')[1])
        else:
            c_ids.append(i)

    cons = Construct.objects.filter(pk__in=c_ids).prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions", "modifications", "deletions", "crystallization__chemical_lists",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")
    proteins = []
    constructs = OrderedDict()
    annotations = {}
    for c in cons:
        # print(c)
        proteins.append(c.protein)
        constructs[c.name] = c.protein.entry_name
        annotations[c.name] = c.schematic()['annotations']

    print(annotations)

    if len(s_ids):
        rs = Residue.objects.filter(protein_conformation__protein__in=proteins, protein_segment__slug__in=s_ids).prefetch_related(
        'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
        'generic_number__scheme', 'display_generic_number__scheme')
    else:
        s_ids = ['N-term','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','ICL4','H8','C-term']
        rs = Residue.objects.filter(protein_conformation__protein__in=proteins).prefetch_related(
        'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
        'generic_number__scheme', 'display_generic_number__scheme')

    print("residues",len(rs))

    distinct_gn = []
    ordered_gn = OrderedDict()
    distinct_segments = []
    overview = OrderedDict()
    segment_length = OrderedDict()
    for s in s_ids:
        overview[s] = OrderedDict()
        segment_length[s] = {'aligned':0, 'before':0,'after':0,'total':0}

    protein_lookup = {}
    print('build stuff')

    segment = ''
    protein = ''

    track_unaligned = {}

    #Find all possible generic numbers, to ensure gaps
    for r in rs.order_by('protein_conformation__id','sequence_number'):
        if segment!=r.protein_segment.slug or protein!=r.protein_conformation.protein.entry_name:
            no_encountered_gn = True
            length = 0
            length_before = 0
            length_after = 0

        segment = r.protein_segment.slug
        protein = r.protein_conformation.protein.entry_name

        if protein not in protein_lookup:
            protein_lookup[protein] = {}
            track_unaligned[protein] = {}

        if segment not in track_unaligned[protein]:
            track_unaligned[protein][segment] = {'before':[],'after':[]}

        if segment not in distinct_segments:
            distinct_segments.append(segment)
            overview[segment] = OrderedDict()

        if r.generic_number:
            no_encountered_gn = False
            gn = r.generic_number.label
            protein_lookup[protein][gn] = {'aa':r.amino_acid,'pos':r.sequence_number}
            gn_sort = gn.split('x')[1]
            gn_sort = float("0."+gn_sort)
            if gn not in distinct_gn:
                distinct_gn.append(gn)
                overview[segment][gn_sort] = [gn,{'aa':'-','pos':''}]
            length += 1
        else:
            if no_encountered_gn:
                track_unaligned[protein][segment]['before'].append({'aa':r.amino_acid,'pos':r.sequence_number})
                length_before += 1
            else:
                track_unaligned[protein][segment]['after'].append({'aa':r.amino_acid,'pos':r.sequence_number})
                length_after += 1

        if len(overview[segment])>segment_length[segment]['aligned']:
            segment_length[segment]['aligned'] = len(overview[segment])
        if length_before>segment_length[segment]['before']:
            segment_length[segment]['before'] = length_before
        if length_after>segment_length[segment]['after']:
            segment_length[segment]['after'] = length_after
        if segment_length[segment]['aligned']+segment_length[segment]['before']+segment_length[segment]['after']>segment_length[segment]['total']:
            segment_length[segment]['total'] = segment_length[segment]['aligned']+segment_length[segment]['before']+segment_length[segment]['after']

    # SORT generic residues to ensure correct order
    gn_list = ""
    ordered_summary = OrderedDict()
    for seg,gns in overview.items():
        ordered_summary[seg] = OrderedDict()
        #GN LIST
        gn_list += """<td class="ali-td ali-residue res-color-X">&nbsp;</td>"""
        if seg!='C-term':
            for _ in range(segment_length[seg]['before']):
                gn_list += """<td class="ali-td">&nbsp;</td>"""
        for gn in sorted(gns):
            ordered_summary[seg][gns[gn][0]] = {'aa':'-','pos':''}
            gn_list += """<td class="ali-td-generic-num">{}</td>""".format("x"+gns[gn][0].split("x")[1])

        if seg=='C-term':
            for _ in range(segment_length[seg]['before']):
                gn_list += """<td class="ali-td">&nbsp;</td>"""

        for _ in range(segment_length[seg]['after']):
            gn_list += """<td class="ali-td">&nbsp;</td>"""

    alignment = OrderedDict()
    alignment_print_sequence = ""
    for c,p in constructs.items():
        alignment[c] = copy.deepcopy(ordered_summary)
        alignment_print_sequence += '<tr>'
        for seg,gns in alignment[c].items():

            if p not in track_unaligned:
                track_unaligned[p] = {seg: {'before':[],'after':[]}}

            if p not in protein_lookup:
                protein_lookup[p] = {}

            if seg not in track_unaligned[p]:
                track_unaligned[p][seg] = {'before':[],'after':[]}

            alignment_print_sequence += """<td class="ali-td ali-residue res-color-_">&nbsp;</td>"""

            if seg!='C-term':
                for _ in range(segment_length[seg]['before']-len(track_unaligned[p][seg]['before'])):
                    alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                    -</td>"""

            for aa in track_unaligned[p][seg]['before']:
                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''

                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],annotation_text,aa['aa'])
            for gn, aa in gns.items():
                if gn in protein_lookup[p]:
                    aa = protein_lookup[p][gn]
                    alignment[c][seg][gn] = aa

                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}<br>SCHEME: {}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],gn,annotation_text,aa['aa'])

            for aa in track_unaligned[p][seg]['after']:
                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],annotation_text,aa['aa'])

            for _ in range(segment_length[seg]['after']-len(track_unaligned[p][seg]['after'])):
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                -</td>"""
            if seg=='C-term':
                for _ in range(segment_length[seg]['before']-len(track_unaligned[p][seg]['before'])):
                    alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                    -</td>"""

        alignment_print_sequence += '</tr>'

    print('done',len(alignment_print_sequence))
    context = {'constructs': constructs,'alignment_print_sequence': alignment_print_sequence, 'segment_length' : segment_length, 'gn_list' : gn_list, 'segments': s_ids, 'c_ids': json.dumps(c_ids)} #, 'alignment_print_sequence': alignment_print_sequence

    return render(request,'align.html',context)
