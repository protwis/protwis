from django.shortcuts import render
from django.conf import settings
from django.core.files import File
from protein.models import ProteinFamily, ProteinAlias, ProteinSet, Protein
from common.views import AbsSettingsSelection
from common.views import AbsSegmentSelection
from common.views import AbsTargetSelection
from common.selection import SelectionItem
import os, shutil, subprocess
import uuid

from phylogenetic_trees.PrepareTree import *
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
            'label': 'Continue to next step',
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

class ButtonSelection():
    """
    """

def render_tree(request):
    Additional_info = {}
    Tree = PrepareTree()
    sets = ProteinSet.objects.all()
    crysts=[]
    for n in sets:
        if n.id==1:
            for prot in n.proteins.all():
                crysts.append(prot.accession)
    # get the user selection from session
    a=Alignment()
    simple_selection=request.session.get('selection', False)
    cons_prots = []
################################## FOR BUILDING STATISTICS ONLY
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
    if len(a.proteins) == 0:
        return render(request, 'phylogenetic_trees/alignment.html', {'phylo': 'None', 'branch':None, 'ttype': None, 'count':None, 'leg':None, 'b':None, 'default':None })  
    a.load_segments_from_selection(simple_selection)
    bootstrap,UPGMA,branches,ttype = map(int,simple_selection.tree_settings)
    bootstrap=10^int(bootstrap)
    if bootstrap==1:
        bootstrap=0 
    #### Create an alignment object
    a.build_alignment()
    a.calculate_statistics()
    a.calculate_similarity()
    total = len(a.proteins)
    lengths= list(map(len,a.proteins[0].alignment)) 
    total_length = len(lengths)-1+sum(lengths)
    families = ProteinFamily.objects.all()
    famdict = {}
    for n in families:
        famdict[Tree.trans_0_2_A(n.slug)]=n.name
   
    dirname = unique_filename = uuid.uuid4()
    os.mkdir('/tmp/%s' %dirname)
    infile = open('/tmp/%s/infile' %dirname,'w')
    infile.write('    '+str(total)+'    '+str(total_length)+'\n')
    family = {}
    Additional_info['crystal']={'proteins':[],'colours':['#6dcde1','#6dcde1']}
    ####Get additional protein information
    for n in a.proteins:
        link = n.protein.entry_name
        name = n.protein.name.replace('<sub>','').replace('</sub>','').replace('<i>','').replace('</i>','')
        if '&' in name and ';' in name:
            name = name.replace('&','').replace(';',' ')
        fam = n.protein.family.slug
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
        fam = Tree.trans_0_2_A(fam)
        if acc in crysts:
            Additional_info['crystal']['proteins'].append(fam)
        if len(name)>25:
            name=name[:25]+'...'
        family[acc] = {'name':name,'family':fam,'description':desc,'species':spec,'class':'','ligand':'','type':'','link': link}
        ####Write PHYLIP input
        sequence = ''
        for chain in n.alignment:
            for residue in n.alignment[chain]:
                sequence += residue[2].replace('_','-')
            sequence += '-'
        infile.write(acc+' '*9+sequence+'\n')
    infile.close()
    ####Run bootstrap
    bootstrap = 0 #TODO

    if bootstrap:
    ### Write phylip input options
        inp = open('/tmp/%s/temp' %dirname,'w')
        inp.write('\n'.join(['r',str(bootstrap),'y','77','y'])+'\n')
        inp.close()
    ###
        subprocess.check_output(['phylip seqboot<temp'], shell=True, cwd = '/tmp/%s' %dirname)

        os.rename('/tmp/%s/infile' %dirname, 'sequence')
        os.rename('/tmp/%s/outfile' %dirname, 'infile')
    ### Write phylip input options
    
    inp = open('/tmp/%s/temp' %dirname,'w')
    if bootstrap:
        inp.write('\n'.join(['m','d',str(bootstrap),'y'])+'\n')
    else:
        inp.write('y\n')
    inp.close()
    ###
    subprocess.check_output(['phylip protdist<temp>>log'], shell=True, cwd = '/tmp/%s' %dirname)
    os.rename('/tmp/%s/infile' %dirname, '/tmp/%s/dupa' %dirname)
    os.rename('/tmp/%s/outfile' %dirname, '/tmp/%s/infile' %dirname)
    inp = open('/tmp/%s/temp' %dirname,'w')
    if bootstrap:
    ### Write phylip input options
        if UPGMA:
            inp.write('\n'.join(['N','m',str(bootstrap),'111','y'])+'\n')
        else:
            inp.write('\n'.join(['m',str(bootstrap),'111','y'])+'\n')
    else:
        if UPGMA:
            inp.write('N\ny\n')
        else:
            inp.write('y\n')
    inp.close()
    ### 
    subprocess.check_output(['phylip neighbor<temp'], shell=True, cwd = '/tmp/%s' %dirname)
    os.remove('/tmp/%s/infile' %dirname)
    os.rename('/tmp/%s/outfile' %dirname, '/tmp/%s/infile' %dirname)
    if bootstrap:
        os.rename('/tmp/%s/outtree' %dirname, '/tmp/%s/intree' %dirname)
    ### Write phylip input options
        inp = open('/tmp/%s/temp' %dirname,'w')
        inp.write('y\n')
        inp.close()
        subprocess.check_output(['phylip consense<temp'], shell=True, cwd = '/tmp/%s' %dirname)

    Tree.treeDo('/tmp/%s/outtree' %dirname,branches,family,famdict)
    phylogeny_input = open('/tmp/%s/out.xml' %dirname,'r').read().replace('\n','')
    shutil.rmtree('/tmp/%s' %dirname)
    return render(request, 'phylogenetic_trees/alignment.html', {'phylo': phylogeny_input, 'branch':branches, 'ttype': ttype, 'count':total, 'leg':str(Tree.legend), 'b':Tree.box, 'default':Tree.defaultColours, 'add':Additional_info })
