from django.contrib.postgres.aggregates import ArrayAgg
from django.conf import settings
from django.core.files import File
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt

from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection
from common.views import AbsMiscSelection
from common.selection import SimpleSelection, Selection, SelectionItem
from mutation.models import *
from phylogenetic_trees.PrepareTree import *
from protein.models import ProteinFamily, ProteinAlias, ProteinSet, Protein, ProteinSegment, ProteinGProteinPair

from copy import deepcopy
import json
import math
import os, shutil, subprocess, signal
import uuid

from collections import OrderedDict

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

def kill_phylo(): #FIXME, needs better way of handling this!
    p = subprocess.Popen(['ps', '-A'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    for line in out.splitlines():
        if 'protdist' in str(line):
            pid = int(line.split(None, 1)[0])
            os.kill(pid, signal.SIGKILL)

class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 3
    docs = 'sequences.html#phylogeneric-trees'
    selection_boxes = OrderedDict([
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/phylogenetic_trees/segmentselection',
            'color': 'success',
        },
    }


class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 3
    docs = 'sequences.html#phylogeneric-trees'
    selection_boxes = OrderedDict([
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/phylogenetic_trees/treesettings',
            'color': 'success',
        },
    }

class TreeSettings(AbsMiscSelection):
    step = 3
    number_of_steps = 3
    docs = 'sequences.html#phylogeneric-trees'
    title = 'SELECT TREE OPTIONS'
    description = 'Select options for tree generation in the middle column.\nOnce you have selected your' \
        + ' settings, click the green button.'
    docs = '/documentation/similarities'
    selection_boxes = OrderedDict([
        ('targets', True),
        ('segments', True),
    ])

    buttons = {
        'continue': {
            'label': 'Calculate & draw tree',
            'url': '/phylogenetic_trees/render_v3',
            'color': 'success',
        },
    }
    # buttons = OrderedDict({
    #     'continue_v3': {
    #         'label': 'Draw tree using v3 code',
    #         'url': '/phylogenetic_trees/render_v3',
    #         'color': 'success',
    #     },
    #     'continue_v2': {
    #         'label': 'Draw tree using new code',
    #         'url': '/phylogenetic_trees/render_v2',
    #         'color': 'success',
    #     },
    #     'continue': {
    #         'label': 'Draw tree using previous code',
    #         'url': '/phylogenetic_trees/render',
    #         'color': 'success',
    #     }
    # })
    tree_settings = True


class Treeclass:
    family = {}

    def __init__(self):
        self.Additional_info={"crystal": {"include":"False", "order":6, "colours":{"crystal_true":"#6dcde1","crystal_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":None, "child": None, "name":"Crystals"},
                "class": {"include":"True", "order":0, "colours":{}, "proteins":[], "color_type":"grayscale", "parent":[], "child": ["family,ligand"], "name":"Class"},
                "family": {"include":"True", "order":1, "colours":{}, "proteins":[], "color_type":"spectrum", "parent":[], "child": ["ligand"], "name":"Ligand type"},
                "ligand": {"include":"True", "order":2, "colours":{}, "proteins":[], "color_type":"spectrum", "parent":["family","class"], "child": [], "name":"Receptor type"},
                "mutant": {"include":"False", "order":3, "colours":{"mutant_true":"#6dcde1","mutant_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":[], "child": ["mutant_plus","mutant_minus"], "name":"Mutated proteins"},
                "mutant_plus": {"include":"False", "order":4, "colours":{"mutant_plus_true":"#6dcde1","mutant_plus_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":"mutant", "child": [], "name":"Positive affinity mutants"},
                "mutant_minus": {"include":"False", "order":5, "colours":{"mutant_minus_true":"#6dcde1","mutant_minus_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":"mutant", "child": [], "name":"Negative affinity mutants"}
                }
        self.buttons = [(x[1]['order'],x[1]['name']) for x in sorted(self.Additional_info.items(), key= lambda x: x[1]['order']) if x[1]['include']=='True']
        self.family = {}
        self.phylip = None
        self.outtree = None
        self.dir = ''

    def Prepare_file(self, request,build=False):
        self.Tree = PrepareTree(build)
        a=Alignment()

        sets = ProteinSet.objects.all()
        #### Get additional data ####
        crysts=[]
        for n in sets:
            if n.id==1:
                for prot in n.proteins.all():
                    crysts.append(prot.entry_name)


        #############################
        # get the user selection from session
        if build != False:
        ################################## FOR BUILDING STATISTICS ONLY##########################
            build_proteins=[]
            if build == '001':
                cons_prots = []
                for prot in Protein.objects.filter(sequence_type__slug='consensus', species_id=1):
                    if prot.family.slug.startswith('001') and len(prot.family.slug.split('_'))==3:
                        build_proteins.append(prot)
                for set in sets:
                    if set.id==1:
                        for prot in set.proteins.all():
                            if prot.family.slug.startswith('001_') and prot.species.latin_name=='Homo sapiens':
                                build_proteins.append(prot)
            else:
                for prot in Protein.objects.filter(sequence_type__slug='wt', species_id=1):
                    if prot.family.slug.startswith(build):
                        build_proteins.append(prot)
            a.load_proteins(build_proteins)
            segments = ProteinSegment.objects.all()
            a.load_segments(segments)
            self.bootstrap,self.UPGMA,self.branches,self.ttype=[0,1,1,0]
        ##################################################################
        else:
            simple_selection=request.session.get('selection', False)
            a.load_proteins_from_selection(simple_selection)
            a.load_segments_from_selection(simple_selection)
            self.bootstrap,self.UPGMA,self.branches,self.ttype = map(int,simple_selection.tree_settings)
        if self.bootstrap!=0:
            self.bootstrap=pow(10,self.bootstrap)
        #### Create an alignment object
        a.build_alignment()
        a.calculate_statistics()
        a.calculate_similarity()
        self.total = len(a.proteins)
        total_length = 0
        for chain in a.proteins[0].alignment:
            total_length += len(a.proteins[0].alignment[chain])
        families = ProteinFamily.objects.all()
        self.famdict = {}
        for n in families:
            self.famdict[self.Tree.trans_0_2_A(n.slug)]=n.name
        dirname = unique_filename = uuid.uuid4()
        os.mkdir('/tmp/%s' %dirname)
        infile = open('/tmp/%s/infile' %dirname,'w')
        infile.write('    '+str(self.total)+'    '+str(total_length)+'\n')
        if len(a.proteins) < 3:
            return 'More_prots',None, None, None, None,None,None,None,None
        ####Get additional protein information
        accesions = {}
        for n in a.proteins:
            fam = self.Tree.trans_0_2_A(n.protein.family.slug)
            if n.protein.sequence_type.slug == 'consensus':
                fam+='_CON'
            entry_name = n.protein.entry_name
            name = n.protein.name.replace('<sub>','').replace('</sub>','').replace('<i>','').replace('</i>','')
            if '&' in name and ';' in name:
                name = name.replace('&','').replace(';',' ')
            acc = n.protein.accession
            if acc:
                acc = acc.replace('-','_')
            else:
                acc = link.replace('-','_')[:6]
            spec = str(n.protein.species)
            fam += '_'+n.protein.species.common_name.replace(' ','_').upper()
            desc = name
            if entry_name in crysts:
                if not fam in self.Additional_info['crystal']['proteins']:
                    self.Additional_info['crystal']['proteins'].append(fam)
            if len(name)>25:
                name=name[:25]+'...'
            self.family[entry_name] = {'name':name,'family':fam,'description':desc,'species':spec,'class':'','accession':acc,'ligand':'','type':'','link': entry_name}
            accesions[acc]=entry_name
            ####Write PHYLIP input
            sequence = ''
            for chain in n.alignment:
                for residue in n.alignment[chain]:
                    sequence += residue[2].replace('_','-')
            infile.write(acc+' '*9+sequence+'\n')
        infile.close()

        ####Run bootstrap
        if self.bootstrap:
        ### Write phylip input options
            inp = open('/tmp/%s/temp' %dirname,'w')
            inp.write('\n'.join(['r',str(self.bootstrap),'y','77','y'])+'\n')
            inp.close()
        ###
            try:
                subprocess.check_output(['phylip seqboot<temp'], shell=True, cwd = '/tmp/%s' %dirname, timeout=300)
                os.rename('/tmp/%s/outfile' %dirname, '/tmp/%s/infile' %dirname)
            except:
                kill_phylo() #FIXME, needs better way of handling this!
                return "too big","too big","too big","too big","too big","too big","too big","too big","too big"

        ### Write phylip input options
        inp = open('/tmp/%s/temp' %dirname,'w')
        if self.bootstrap:
            inp.write('\n'.join(['m','d',str(self.bootstrap),'y'])+'\n')
        else:
            inp.write('y\n')
        inp.close()
        ###
        try:
            subprocess.check_output(['phylip protdist<temp>>log'], shell=True, cwd = '/tmp/%s' %dirname, timeout=300)
        except:
            kill_phylo() #FIXME, needs better way of handling this!
            return "too big","too big","too big","too big","too big","too big","too big","too big","too big"
        os.rename('/tmp/%s/infile' %dirname, '/tmp/%s/dupa' %dirname)
        os.rename('/tmp/%s/outfile' %dirname, '/tmp/%s/infile' %dirname)
        inp = open('/tmp/%s/temp' %dirname,'w')
        if self.bootstrap:
        ### Write phylip input options
            if self.UPGMA:
                inp.write('\n'.join(['N','m',str(self.bootstrap),'111','y'])+'\n')
            else:
                inp.write('\n'.join(['m',str(self.bootstrap),'111','y'])+'\n')
        else:
            if self.UPGMA:
                inp.write('N\ny\n')
            else:
                inp.write('y\n')
        inp.close()
        ###
        try:
            subprocess.check_output(['phylip neighbor<temp'], shell=True, cwd = '/tmp/%s' %dirname, timeout=300)
        except:
            kill_phylo() #FIXME, needs better way of handling this!
            return "too big","too big","too big","too big","too big","too big","too big","too big"
        if self.bootstrap:
            os.rename('/tmp/%s/outfile' %dirname, '/tmp/%s/infile' %dirname)
            os.rename('/tmp/%s/outtree' %dirname, '/tmp/%s/intree' %dirname)
           ### Write phylip input options
            inp = open('/tmp/%s/temp' %dirname,'w')
            inp.write('y\n')
            inp.close()
        ###
            try:
                subprocess.check_output(['phylip consense<temp'], shell=True, cwd = '/tmp/%s' %dirname, timeout=300)
            except:
                kill_phylo() #FIXME, needs better way of handling this!
                return "too big","too big","too big","too big","too big","too big","too big","too big"
        self.phylip = open('/tmp/%s/outtree' %dirname).read()
        for acc in accesions.keys():
            self.phylip=self.phylip.replace(acc,accesions[acc])
#        self.phylogeny_output = self.phylip
        self.outtree = open('/tmp/%s/outfile' %dirname).read().lstrip()
        phylogeny_input = self.get_phylogeny('/tmp/%s/' %dirname)
        shutil.rmtree('/tmp/%s' %dirname)

        if build != False:
            open('static/home/images/'+build+'_legend.svg','w').write(str(self.Tree.legend))
            open('static/home/images/'+build+'_tree.xml','w').write(phylogeny_input)
        else:
            return phylogeny_input, self.branches, self.ttype, self.total, str(self.Tree.legend), self.Tree.box, self.Additional_info, self.buttons, a.proteins

    def get_phylogeny(self, dirname):

        self.Tree.treeDo(dirname, self.phylip,self.branches,self.family,self.Additional_info, self.famdict)
        phylogeny_input = open('%s/out.xml' %dirname,'r').read().replace('\n','')
        return phylogeny_input

    def get_data(self):
        return self.branches, self.ttype, self.total, str(self.Tree.legend), self.Tree.box, self.Additional_info, self.buttons


def get_buttons(request):
    Tree_class=request.session['Tree']
    buttons = [(x[1]['order'],x[1]['name']) for x in sorted(Tree_class.Additional_info.items(), key= lambda x: x[1]['order']) if x[1]['include']=='True']
    return render(request, 'phylogenetic_trees/ring_buttons.html', {'but':buttons })

def modify_tree(request):
    try:
        shutil.rmtree('/tmp/modify')
    except:
        pass
    arg = request.GET.getlist('arg[]')
    value = request.GET.getlist('value[]')
    Tree_class=request.session['Tree']
    for n in range(len(arg)):
        Tree_class.Additional_info[arg[n].replace('_btn','')]['include']=value[n]
    request.session['Tree']=Tree_class
    os.mkdir('/tmp/modify')
    phylogeny_input = Tree_class.get_phylogeny('/tmp/modify')
    branches, ttype, total, legend, box, Additional_info, buttons=Tree_class.get_data()
    shutil.rmtree('/tmp/modify')
    if ttype == '1':
        float(total)/4*100
    else:
        count = 1900 - 1400/math.sqrt(float(total))

    return render(request, 'phylogenetic_trees/main.html', {'phylo': phylogeny_input, 'branch':branches, 'ttype': ttype, 'count':count, 'leg':legend, 'b':box, 'add':Additional_info, 'but':buttons, 'phylip':Tree_class.phylip, 'outtree':Tree_class.outtree})

# DEPRECATED CODE - can be cleaned up
def render_tree(request):
    Tree_class=Treeclass()
    phylogeny_input, branches, ttype, total, legend, box, Additional_info, buttons, proteins=Tree_class.Prepare_file(request)
    if phylogeny_input == 'too big':
        return render(request, 'phylogenetic_trees/too_big.html')

    if phylogeny_input == 'More_prots':
        return render(request, 'phylogenetic_trees/warning.html')

    if ttype == '1':
        float(total)/4*100
    else:
        count = 1900 - 1400/math.sqrt(float(total))

    request.session['Tree']=Tree_class
    return render(request, 'phylogenetic_trees/alignment.html', {'phylo': phylogeny_input, 'branch':branches, 'ttype': ttype, 'count':count, 'leg':legend, 'b':box, 'add':Additional_info, 'but':buttons, 'phylip':Tree_class.phylip, 'outtree':Tree_class.outtree })

# DEPRECATED CODE - can be cleaned up
def render_tree_v2(request):
    Tree_class=Treeclass()
    phylogeny_input, branches, ttype, total, legend, box, Additional_info, buttons, proteins=Tree_class.Prepare_file(request)
    if phylogeny_input == 'too big':
        return render(request, 'phylogenetic_trees/too_big.html')

    if phylogeny_input == 'More_prots':
        return render(request, 'phylogenetic_trees/warning.html')

    if ttype == '1':
        float(total)/4*100
    else:
        count = 1900 - 1400/math.sqrt(float(total))


    protein_data = []

    #FIXME remove
    import random
    for pc in proteins:
        v = {}
        p = pc.protein
        v['name'] = p.entry_name
        v['GPCR_class'] = p.family.parent.parent.parent.name
        v['selectivity'] = ["Gq/G11 family"]
        v['ligand_type'] = p.family.parent.parent.name
        v['coverage'] = random.uniform(0, 1)
        v['receptor_page'] = ''
        print(v)
        protein_data.append(v)


    request.session['Tree']=Tree_class
    context = {}
    context['phylip'] = Tree_class.phylip.replace('\n', '')
    context['protein_data'] = protein_data
    return render(request, 'phylogenetic_trees/display.html', context)

@csrf_exempt
def signature_selection(request):
    # create full selection and import simple selection (if it exists)
    simple_selection = request.session.get('selection', False)
    selection_pos = Selection()
    selection_pos.importer(deepcopy(simple_selection))
    selection_pos.clear('targets')

    selection_neg = Selection()
    selection_neg.importer(deepcopy(simple_selection))
    selection_neg.clear('targets')

    if 'group1' in request.POST and 'group2' in request.POST:
        up_names = request.POST['group1'].split('\r')
        for up_name in up_names:
            try:
                selection_object = SelectionItem('protein', Protein.objects.get(entry_name=up_name.strip().lower()))
                selection_pos.add('targets', 'protein', selection_object)
            except:
                continue

        down_names = request.POST['group2'].split('\r')
        for down_name in down_names:
            try:
                selection_object = SelectionItem('protein', Protein.objects.get(entry_name=down_name.strip().lower()))
                selection_neg.add('targets', 'protein', selection_object)
            except:
                continue

        # Set both the positive and negative target selections
        request.session['targets_pos'] = selection_pos.exporter()

        request.session['selection'] = selection_neg.exporter()

        return JsonResponse({"response": "ok"})
    else:
        return JsonResponse({"response": "error"})

def render_tree_v3(request):
    Tree_class=Treeclass()

    phylogeny_input, branches, ttype, total, legend, box, Additional_info, buttons, proteins=Tree_class.Prepare_file(request)
    if phylogeny_input == 'too big':
        return render(request, 'phylogenetic_trees/too_big.html')

    if phylogeny_input == 'More_prots':
        return render(request, 'phylogenetic_trees/warning.html')

    # if ttype == '1':
    #     float(total)/4*100
    # else:
    #     count = 1900 - 1400/math.sqrt(float(total))

    #protein_data = []
    #
    #FIXME remove
    # import random
    # for pc in proteins:
    #     v = {}
    #     p = pc.protein
    #     v['name'] = p.entry_name
    #     v['GPCR_class'] = p.family.parent.parent.parent.name
    #     v['selectivity'] = ["Gq/G11 family"]
    #     v['ligand_type'] = p.family.parent.parent.name
    #     v['coverage'] = random.uniform(0, 1)
    #     v['receptor_page'] = ''
    #     print(v)
    #     protein_data.append(v)

    request.session['Tree'] = Tree_class

    # output dictionary
    data = {}
    data['tree'] = Tree_class.phylip.replace('\n', '')
    # context['protein_data'] = protein_data


    protein_entries = []
    for pc in proteins:
        protein_entries.append(pc.protein.entry_name)

    # load all
    cluster_method = 0

    # Collect structure annotations
    protein_annotations = {}

    # Grab all annotations and all the ligand role when present in aggregates
    # NOTE: we can probably remove the parent step and go directly via family in the query
    annotations = Protein.objects.filter(entry_name__in=protein_entries) \
                    .values_list('entry_name', 'name', 'family__parent__name', 'family__parent__parent__name', 'family__parent__parent__parent__name', 'family__slug', 'name')

    protein_slugs = set()
    for an in annotations:
        protein_annotations[an[0]] = list(an[1:])

        # add slug to lists
        slug = protein_annotations[an[0]][4]
        protein_slugs.add(slug)

    data['annotations'] = protein_annotations

    # Grab G-protein coupling profile for all receptors covered by the selection
    # TODO: make general cache(s) for these kinds of data
    selectivitydata = {}
    coupling = ProteinGProteinPair.objects.filter(protein__family__slug__in=protein_slugs, source="GuideToPharma").values_list('protein__family__slug', 'transduction').annotate(arr=ArrayAgg('g_protein__name'))

    for pairing in coupling:
        if pairing[0] not in selectivitydata:
            selectivitydata[pairing[0]] = {}
        selectivitydata[pairing[0]][pairing[1]] = pairing[2]

    data['Gprot_coupling'] = selectivitydata

    context = {}
    context["data"] = json.dumps(data)

    return render(request, 'phylogenetic_trees/phylo_tree.html', context)
