from django.db.models import Count, Avg, Min, Max
from collections import defaultdict
from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import TemplateView, View
from protein.models import Protein, Species, ProteinFamily
from ligand.models import Ligand, AssayExperiment, LigandProperities 
# Create your views here.

    
def LigandBrowser(request):
    
    ligands = AssayExperiment.objects.values('protein__entry_name', 'protein__family__name', 'protein__family__parent__name', 'protein__family__parent__parent__name', 'protein__family__parent__parent__parent__name', 'protein__family__slug', 'protein__family__parent__slug', 'protein__family__parent__parent__slug', 'protein__family__parent__parent__parent__slug', 'protein__species__common_name').annotate(num_ligands=Count('ligand', distinct=True))
    
    print(ligands)
    
    #od = collections.OrderedDict(sorted(ligands.items()))
    #
    #context = {'ligands': od}
    context = {'ligands': ligands}

    return render(request, 'assay/ligand_browser.html', context)

 



def l_detail(request, ligand__id):
    print (ligand__id)
    ls = AssayExperiment.objects.filter(ligand__properities__web_links__index=ligand__id)
    ligands = ls.values(
        'protein__entry_name', 'protein__family__name', 'protein__family__parent__name',
        'protein__family__parent__parent__name', 'protein__family__parent__parent__parent__name',
        'assay_type','standard_type', 'standard_units').annotate(num_records = Count('protein__entry_name')).order_by('protein__entry_name')
    
    
    ps =ls.values('standard_type', 'standard_value', 'assay_type', 'protein__entry_name','published_value', 'published_type', 'assay_id').order_by('protein__entry_name')#.annotate(Avg('standard_value'))#, distinct=True))
    #for p in ps:
    #    #print(p.__dict__)
    #    print (p)#,p['assay_type'], p['standard_type'], p['standard_value'], p['c'])
    ##print ('ps{}'.format(ps))
    #
    p_dict = {}
    p_dict = defaultdict(lambda: defaultdict(lambda:defaultdict(list)))
    for p in ps:
        #print (p)
        p_dict[p['protein__entry_name']][p['assay_type']][p['standard_type']].append(float(p['standard_value']))
    #print (p_dict['5ht1a_human'])
    for i in p_dict:
        for l in p_dict[i]: 
            for d in p_dict[i][l]:
                #print (p_dict[i][l][d])
                #print (len(p_dict[i][l][d]))
                
                av = (sum(p_dict[i][l][d]))/len(p_dict[i][l][d]) #average acros different measuremnts
                
                av = float("{0:.2f}".format(av)) #average acros different measuremnts
                min_val = min(p_dict[i][l][d]) #maximum value
                p_dict[i][l][d]=[min_val,av]
            
    for i in p_dict:
        print(i, p_dict[i])
    #
    av_list = []
    for p in ligands:
        #print (p)
        p['average']=p_dict[p['protein__entry_name']][p['assay_type']][p['standard_type']][1]
        p['min']=p_dict[p['protein__entry_name']][p['assay_type']][p['standard_type']][0]
        if p['protein__entry_name'] not in av_list:
            av_list.append(p)
    print(len(ligands),len(ps))
    
    #for i in av_list:
    #    print (i['protein__entry_name'],i['assay_type'],i['standard_type'],i['average'])    
    ##context = {'ligands': ligands,'assay_b':ls1, 'assay_f':ls2}
    context = {'ligands': av_list, 's':ligand__id}
    
    
    
    return render(request, 'assay/l_detail.html', context)

    
    
def p_details_filtered(request, slug):
    # copy paste ish
    
    # filter
    pass

def p_detail(request, slug):
    print (slug)
    if slug.count('_') == 0 :
        ps = AssayExperiment.objects.filter(protein__family__parent__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
    elif slug.count('_') == 1 and len(slug) == 7:
        ps = AssayExperiment.objects.filter(protein__family__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
    elif slug.count('_') == 2:
        ps = AssayExperiment.objects.filter(protein__family__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
    #elif slug.count('_') == 3:
    elif slug.count('_') == 1 and len(slug) != 7:
        ps = AssayExperiment.objects.filter(protein__entry_name = slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
    
    ps = ps.values('standard_type', 'standard_relation', 'standard_value', 'assay_description', 'assay_type', 'standard_units', 'pchembl_value',
                   'ligand__id', 'ligand__properities_id', 'ligand__properities__web_links__index',
                   'protein__species__common_name', 'protein__entry_name', 'ligand__properities__mw',
                   'ligand__properities__logp', 'ligand__properities__rotatable_bonds', 'ligand__properities__smiles',
                   'ligand__properities__hdon', 'ligand__properities__hacc','protein', 'ligand__properities__ligandvendorlink__vendor__name'
                   ).annotate(num_targets = Count('protein__id', distinct=True))

    if slug.count('_') == 1 and len(slug) == 7:
        f = ProteinFamily.objects.get(slug=slug)      
    else:
        f = slug
    
    
    for p in ps:
        #print (p.standard_value)
        print (p['protein__entry_name'])
    #print(context)
    print ('hi')
    #f = Protein.objects.get(family__slug='slug')http://localhost:8000/ligand/
    #print(f)
    print (len(ps))
    #context = {'proteins': ps, 's':slug}
    context = {'proteins': ps, 's':f}

    #def average(request,):
    #    p_dict = defaultdict(lambda: defaultdict(lambda:defaultdict(list)))
    #    p_dict[p['ligand__properities__web_links__index']][p['assay_type']][p['standard_type']].append(float(p['standard_value']))
    #    for i in p_dict:                                                   
    #        for d in p_dict[i][l]:
    #            av = sum(p_dict[i][l][d])/len(p_dict[i][l][d])
    #            p_dict[i][l][d].append(av)
    #    
    #            
        
    return render(request, 'assay/p_detail.html', context)



#ps = AssayExperiment.objects.filter(protein__entry_name = slug, ligand__properities__web_links__web_resource_id = 10).values('standard_type', 'standard_relation', 'standard_value', 'assay', 'standard_units','ligand__id', 'ligand__properities__web_links__index','protein__species__common_name', 'protein__entry_name').annotate(num_targets = Count('protein'))#, distinct=True))
    #ligands=LigandProperities.objects.filter(ligand__id = ps.ligand__id)
    #for record in ps:
    #    ls.append = LigandProperities.objects.filter(properities_id=record['ligand__id'])


    #ps = AssayExperiment.objects.filter(protein__entry_name=slug, ligand__properities__web_links__web_resource_id = 10).values('standard_type', 'standard_relation', 'standard_value', 'assay', 'standard_units','ligand__id', 'ligand__properities__web_links__index','protein__species__common_name').annotate(num_targets = Count('protein'))#, distinct=True))
    #ps = AssayExperiment.objects.filter(protein__entry_name = slug).values('lignds').annotate(num_ligands=Count('ligand', distinct = True))
    #ps = AssayExperiment.objects.filter(protein__entry_name='5ht1e_human').values('ligand__properities__mw', 'ligand__properities__rotatable_bonds', 'standard_type', 'standard_relation', 'standard_value', 'assay')
    





 #protein_links = p.web_links.all().distinct('web_resource__slug')   'ligand__properities__web_links__index'
  #ps = AssayExperiment.objects.filter(protein__entry_name='5ht1e_human').values('ligand__properities__mw', 'ligand__properities__rotatable_bonds', 'standard_type', 'standard_relation', 'standard_value', 'assay''ligand__properities__web_links__index')


#def ligand_browser (request):
#    pls = models.AssayExperiment.objects.values('protein__entry_name','protein__family__name').annotate(num_ligands=Count('ligand', distinct=True))
##    pls = models.AssayExperiment.objects.values('protein__entry_name','protein__family__name', 'protein__family__parent__name', 'protein__family__parent__parent__name').annotate(num_ligands=Count('ligand', distinct=True))
##    for pl in pls:
##        print(pl)
##    m.objects.values('p').annotate(number_of_days=Count('date', distinct=True))
##    p = Protein.objects.prefetch_related('web_links__web_resource').filter(sequence_type__slug='wt')
#    return render(request,'assay/ligand_browser.html', pls) # ot pls