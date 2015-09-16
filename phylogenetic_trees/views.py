from django.shortcuts import render
from django.conf import settings
from django.core.files import File
from protein.models import ProteinFamily, ProteinAlias, ProteinSet, Protein
from common.views import AbsSettingsSelection
from common.views import AbsSegmentSelection
from common.views import AbsTargetSelection
from common.selection import SelectionItem
from mutation.models import *

import os, shutil, subprocess
import uuid
from phylogenetic_trees.PrepareTree2 import *

# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from collections import OrderedDict
#TODO
class TreeSettings(AbsSettingsSelection):
    step = 3
    number_of_steps = 3
    docs = '/docs/phylogenetic_trees'
    selection_boxes = OrderedDict([
        ('tree_settings', False),
        ('segments', True),
        ('targets', True),
    ])
    buttons = {
        'continue': {
            'label': 'Draw tree',
            'url': '/phylogenetic_trees/render',
            'color': 'success',
        },
   }

class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 3
    docs = '/docs/phylogenetic_trees'
    selection_boxes = OrderedDict([
        ('tree_settings', False),
        ('segments', True),
        ('targets', True),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/phylogenetic_trees/treesettings',
            'color': 'success',
        },
    }


class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 3
    docs = '/docs/phylogenetic_trees'
    selection_boxes = OrderedDict([
        ('tree_settings', False),
        ('segments', False),
        ('targets', True),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/phylogenetic_trees/segmentselection',
            'color': 'success',
        },
    }

class Treeclass:
    family = {}
    phylip = None
    dir = ''
    Additional_info={"crystal": {"include":"False", "order":6, "colours":{"crystal_true":"#6dcde1","crystal_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":None, "child": None, "name":"Crystals"},
            "class": {"include":"True", "order":0, "colours":{}, "proteins":[], "color_type":"grayscale", "parent":[], "child": ["family,ligand"], "name":"Class"},
            "family": {"include":"True", "order":1, "colours":{}, "proteins":[], "color_type":"spectrum", "parent":[], "child": ["ligand"], "name":"Ligand type"},
            "ligand": {"include":"True", "order":2, "colours":{}, "proteins":[], "color_type":"spectrum", "parent":["family","class"], "child": [], "name":"Ligand name"},
            "mutant": {"include":"False", "order":3, "colours":{"mutant_true":"#6dcde1","mutant_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":[], "child": ["mutant_plus","mutant_minus"], "name":"Mutated proteins"},
            "mutant_plus": {"include":"False", "order":4, "colours":{"mutant_plus_true":"#6dcde1","mutant_plus_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":"mutant", "child": [], "name":"Positive affinity mutants"},
            "mutant_minus": {"include":"False", "order":5, "colours":{"mutant_minus_true":"#6dcde1","mutant_minus_false":"#EEE"}, "color_type":"single", "proteins":[], "parent":"mutant", "child": [], "name":"Negative affinity mutants"}
            }
    buttons = [(x[1]['order'],x[1]['name']) for x in sorted(Additional_info.items(), key= lambda x: x[1]['order']) if x[1]['include']=='True']

    def Prepare_file(self, request):
        self.Tree = PrepareTree()
        sets = ProteinSet.objects.all()
        #### Get additional data ####
        crysts=[]
        for n in sets:
            if n.id==1:
                for prot in n.proteins.all():
                    crysts.append(prot.accession)

            
        #############################
        # get the user selection from session
        a=Alignment()
        simple_selection=request.session.get('selection', False)
        ################################## FOR BUILDING STATISTICS ONLY
        #cons_prots = []
        #    for n in Protein.objects.filter(sequence_type_id=3):
        #        if n.family.slug.startswith('001') and len(n.family.slug.split('_'))==2:
        #            cons_prots.append(n)
        #    for n in simple_selection.targets:
        #        if n.item.family.slug.startswith('001_'):
        #            continue
        #        else:
        #            simple_selection.targets.remove(n)
        #    for n in cons_prots:
        #        simple_selection.targets.append(SelectionItem('protein',n))
        #####################################################
        a.load_proteins_from_selection(simple_selection)
        a.load_segments_from_selection(simple_selection)
        self.bootstrap,self.UPGMA,self.branches,self.ttype = map(int,simple_selection.tree_settings)
        self.bootstrap=10^int(self.bootstrap)
        if self.bootstrap==1:
            self.bootstrap=0 
        #### Create an alignment object
        a.build_alignment()
        a.calculate_statistics()
        a.calculate_similarity()
        self.total = len(a.proteins)
        lengths= list(map(len,a.proteins[0].alignment)) 
        total_length = len(lengths)-1+sum(lengths)
        families = ProteinFamily.objects.all()
        self.famdict = {}
        for n in families:
            self.famdict[self.Tree.trans_0_2_A(n.slug)]=n.name
        dirname = unique_filename = uuid.uuid4()
        os.mkdir('/tmp/%s' %dirname)
        infile = open('/tmp/%s/infile' %dirname,'w')
        infile.write('    '+str(self.total)+'    '+str(total_length)+'\n')
        ####Get additional protein information
        muts = {}
        
        for n in a.proteins:
            fam = self.Tree.trans_0_2_A(n.protein.family.slug)

            #### Mutations ######
            #muts_db = MutationExperiment.objects.filter(protein_id=n.id)
            #for a in muts_db:
            #    print(a.protein,a.residue,a.mutation,a.ligand)
            #    wt = int(a.wt_value)
            #    mu = int(a.mu_value)
            #    if wt < mu:
            #        muts[fam]='negative'
            #    if wt > mu:
            #        muts[fam]='positive'
            #    if wt == mu:
            #        muts[fam]='neutral'
            #exit()

            link = n.protein.entry_name
            name = n.protein.name.replace('<sub>','').replace('</sub>','').replace('<i>','').replace('</i>','')
            if '&' in name and ';' in name:
                name = name.replace('&','').replace(';',' ')
            acc = n.protein.accession
            if acc:
                acc = acc.replace('-','_')
            else:
                acc = link.replace('-','_')[:6]
            spec = str(n.protein.species)
            try:
                desc = str(ProteinAlias.objects.filter(protein__in=[n.id])[0])
            except IndexError:
                desc = ''
            if acc in crysts:
                if not fam in self.Additional_info['crystal']['proteins']:
                    self.Additional_info['crystal']['proteins'].append(fam)
            if len(name)>25:
                name=name[:25]+'...'
            self.family[acc] = {'name':name,'family':fam,'description':desc,'species':spec,'class':'','ligand':'','type':'','link': link}
            ####Write PHYLIP input
            sequence = ''
            for chain in n.alignment:
                for residue in n.alignment[chain]:
                    sequence += residue[2].replace('_','-')
                sequence += '-'
            infile.write(acc+' '*9+sequence+'\n')
        infile.close()
        ####Run bootstrap
        if self.bootstrap:
        ### Write phylip input options
            inp = open('/tmp/%s/temp' %dirname,'w')
            inp.write('\n'.join(['r',str(self.bootstrap),'y','77','y'])+'\n')
            inp.close()
        ###
            subprocess.check_output(['phylip seqboot<temp'], shell=True, cwd = '/tmp/%s' %dirname)
            os.rename('/tmp/%s/outfile' %dirname, '/tmp/%s/infile' %dirname)
        ### Write phylip input options
        inp = open('/tmp/%s/temp' %dirname,'w')
        if self.bootstrap:
            inp.write('\n'.join(['m','d',str(self.bootstrap),'y'])+'\n')
        else:
            inp.write('y\n')
        inp.close()
        ###
        subprocess.check_output(['phylip protdist<temp>>log'], shell=True, cwd = '/tmp/%s' %dirname)
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
        subprocess.check_output(['phylip neighbor<temp'], shell=True, cwd = '/tmp/%s' %dirname)
        os.remove('/tmp/%s/infile' %dirname)
        os.rename('/tmp/%s/outfile' %dirname, '/tmp/%s/infile' %dirname)
        if self.bootstrap:
           ### Write phylip input options
            os.rename('/tmp/%s/outtree' %dirname, '/tmp/%s/intree' %dirname)
            inp = open('/tmp/%s/temp' %dirname,'w')
            inp.write('y\n')
            inp.close()
        ###
            subprocess.check_output(['phylip consense<temp'], shell=True, cwd = '/tmp/%s' %dirname)
        self.phylip = open('/tmp/%s/outtree' %dirname).read()
        phylogeny_input = self.get_phylogeny('/tmp/%s/' %dirname)
        #shutil.rmtree('/tmp/%s' %dirname)
        return phylogeny_input, self.branches, self.ttype, self.total, str(self.Tree.legend), self.Tree.box, self.Additional_info, self.buttons
        
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
    arg = request.GET.getlist('arg[]')
    value = request.GET.getlist('value[]')
    Tree_class=request.session['Tree']
    for n in range(len(arg)):
        Tree_class.Additional_info[arg[n].replace('_btn','')]['include']=value[n]
    request.session['Tree']=Tree_class
    phylogeny_input = Tree_class.get_phylogeny('/tmp')
    branches, ttype, total, legend, box, Additional_info, buttons=Tree_class.get_data()
    
    return render(request, 'phylogenetic_trees/main.html', {'phylo': phylogeny_input, 'branch':branches, 'ttype': ttype, 'count':total, 'leg':legend, 'b':box, 'add':Additional_info, 'but':buttons })

def render_tree(request):
    Tree_class=Treeclass()
    phylogeny_input, branches, ttype, total, legend, box, Additional_info, buttons=Tree_class.Prepare_file(request)
    request.session['Tree']=Tree_class
    return render(request, 'phylogenetic_trees/alignment.html', {'phylo': phylogeny_input, 'branch':branches, 'ttype': ttype, 'count':total, 'leg':legend, 'b':box, 'add':Additional_info, 'but':buttons })


