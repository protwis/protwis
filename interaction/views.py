from django.shortcuts import render
from django.http import HttpResponse
#from interaction.functions import *
# Create your views here.
from subprocess import call
from .forms import PDBform
from django import forms
from django.core.servers.basehttp import FileWrapper
from django.db.models import Count, Min, Sum, Avg

from interaction.models import *
from ligand.models import Ligand
from ligand.models import LigandType
from ligand.models import LigandRole
from structure.models import Structure,PdbData,Rotamer,Fragment
from protein.models import ProteinConformation
from residue.models import Residue
from common.models import WebResource
from common.models import WebLink


from os import path, listdir
from os.path import isfile, join
import yaml
from operator import itemgetter
import datetime
import re

AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

module_dir = path.dirname(__file__)

def regexaa(aa):
    aaPattern = re.compile(r'^(\w{3})(\d+)(\w+)$') 
    result = aaPattern.search(aa).groups() #Splits the string into AA number CHAIN : LEU339A => ('LEU', '339', 'A')
    if result:
        aa = AA[result[0]]
        number = result[1]
        chain = result[2]
        return aa,number,chain
    else:
        return None, None, None

def index(request):

    form = PDBform()

    #context = {}
    return render(request,'interaction/index.html',{'form': form})

def list(request):
    form = PDBform()
    structures = StructureLigandInteraction.objects.distinct('structure').all()
    print(structures)
    for structure in structures:
        print(structure.structure.pdb_code.index)
    #objects = Model.objects.filter(id__in=object_ids)
    #context = {}
    return render(request,'interaction/list.html',{'form': form, 'structures': structures})

def view(request):
    pdbname = request.GET.get('pdb')
    form = PDBform()
    structures = StructureLigandInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(residuefragmentatom__atomtype='ATOM').annotate(numRes = Count('residuefragmentatom__residuenr', distinct = True)).order_by('-numRes')
    #structures = StructureLigandInteraction.objects.annotate(test = Count(ResidueFragmentAtom))

    return render(request,'interaction/view.html',{'form': form, 'pdbname': pdbname, 'structures': structures})

def ligand(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    form = PDBform()
    fragments = ResidueFragmentInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(ligand__name=ligand).order_by('interaction_type')
    #print(fragments)
    for fragment in fragments:
        print(fragment.fragment)
    #objects = Model.objects.filter(id__in=object_ids)
    #context = {}
    return render(request,'interaction/ligand.html',{'form': form, 'pdbname': pdbname, 'ligand': ligand, 'fragments': fragments})

def calculate(request):     
    if request.method == 'POST':
        form = PDBform(request.POST)
        if form.is_valid():
            pdbname = form.cleaned_data['pdbname']
            call(["python", "interaction/functions.py","-p",pdbname])
            mypath = module_dir+'/temp/results/'+pdbname+'/output'


            results = []
            web_resource, created = WebResource.objects.get_or_create(slug='pdb',url='http://www.rcsb.org/pdb/explore/explore.do?structureId=$index')
            web_link, created = WebLink.objects.get_or_create(web_resource=web_resource,index=pdbname)

            structure=Structure.objects.filter(pdb_code=web_link) 
            if structure.exists():
                structure=Structure.objects.get(pdb_code=web_link)
            else:
                 quit() #quit!

            protein=structure.protein_conformation

            for f in listdir(mypath):
                if isfile(join(mypath,f)):       
                    result = yaml.load(open(mypath+"/"+f, 'rb'))
                    output = result

                    temp = f.replace('.yaml','').split("_")
                    temp.append([output])
                    temp.append(round(output['score'][0][0]))
                    temp.append((output['inchikey']).strip())
                    temp.append((output['smiles']).strip())
                    results.append(temp)

                    ligandtype, created = LigandType.objects.get_or_create(slug="sm", name='Small molecule')

                    ligand, created = Ligand.objects.get_or_create(inchikey=output['inchikey'].strip(),
                        defaults={'name':temp[1], 'smiles':output['smiles'].strip(),'ligand_type':ligandtype})


                    #proteinligand, created = ProteinLigandInteraction.objects.get_or_create(protein=protein,ligand=ligand)

                    
                    ligandrole, created = LigandRole.objects.get_or_create(name='unknown',slug='unknown')

                    f = module_dir+"/temp/results/"+pdbname+"/interaction"+"/"+pdbname+"_"+temp[1]+".pdb"
                    if isfile(f):      
                        pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read()) #does this close the file?
                        print("Found file"+f)
                    else:
                        quit()

                    structureligandinteraction, created = StructureLigandInteraction.objects.get_or_create(ligand=ligand,structure=structure, ligand_role=ligandrole, pdb_file=pdbdata) #, pdb_reference=pdbname <-- max length set to 3? So can't insert ones correctly

                    
                    for interactiontype,interactionlist in output.items():
                        if interactiontype=='hbond' or interactiontype=='hbondplus':
                            for entry in interactionlist:
                                #print(interactiontype,entry)
                                aa = entry[0]
                                aa,pos,chain = regexaa(aa)

                                residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                if residue.exists():
                                    residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                else:
                                    print("Not found residue!",pdbname,pos,aa)
                                    continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                                for pair in entry[3]:
                                    fragment = pair[0] #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                                    f = module_dir+"/temp/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_"+fragment+"_HB.pdb"
                                    if interactiontype=='hbondplus':  f = module_dir+"/temp/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_"+fragment+"_HBC.pdb"
                                    if isfile(f):      
                                        print("Found file"+f)
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
                                        f_in.close();

                                        rotamer_data, created = PdbData.objects.get_or_create(pdb=rotamer_pdb)
                                        rotamer, created = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)

                                        fragment_data, created = PdbData.objects.get_or_create(pdb=fragment_pdb) 
                                        fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=fragment_data)
                                    else:
                                        quit("Could not find "+f)

                                    interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug=interactiontype,name=interactiontype)
                                    fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure=structure,residue=residue,ligand=ligand,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)
                                    
                                    
                        elif interactiontype=='hbond_confirmed':
                            for entry in interactionlist:
                                #print(interactiontype,entry)
                                aa = entry[0]
                                aa,pos,chain = regexaa(aa)
                                interactiontype="HB"+entry[1][0][0]

                                residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                if residue.exists():
                                    residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                else:
                                    print("Not found residue!",pdbname,pos,aa)
                                    continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                                for pair in entry[1]:  
                                    fragment = pair[1] #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                                    #mypath = module_dir+'/temp/results/'+pdbname+'/fragments'
                                    f = module_dir+"/temp/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_"+fragment+"_HB.pdb"
                                    if isfile(f):      
                                        print("Found file"+f)
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
                                        f_in.close();

                                        rotamer_data, created = PdbData.objects.get_or_create(pdb=rotamer_pdb)
                                        rotamer, created = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)

                                        fragment_data, created = PdbData.objects.get_or_create(pdb=fragment_pdb) 
                                        fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=fragment_data)
                                    else:
                                        quit("Could not find "+f)

                                    interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug=interactiontype,name=interactiontype)
                                    fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure=structure,residue=residue,ligand=ligand,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)
                          
                        elif interactiontype=='hydrophobic': #NO FRAGMENT FOR THESE, WHAT TO DO? USE WHOLE LIGAND?
                            for entry in interactionlist:
                                #print(entry)
                                aa = entry[0]
                                aa,pos,chain = regexaa(aa)
                                #fragment = entry[1][0][1] #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                                residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                if residue.exists():
                                    residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                else:
                                    print("Not found residue!",pdbname,pos,aa)
                                    continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                                fragment = '' #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                                #interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug='hydrofob',name=interactiontype)
                                #fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure=structure,residue=residue,ligand=ligand,interaction_type=interaction_type,fragment=fragment)
                        elif interactiontype=='aromaticplus' or interactiontype=='aromatic':
                            for entry in interactionlist:
                                print(entry)
                                aa = entry[0]
                                aa,pos,chain = regexaa(aa)
                                fragment = entry[1]

                                residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                if residue.exists():
                                    residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                else:
                                    print("Not found residue!",pdbname,pos,aa)
                                    continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                                f = module_dir+"/temp/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_aromatic_"+str(entry[1])+".pdb"
                                if isfile(f):      
                                    print("Found file"+f)
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
                                    f_in.close();   

                                    rotamer_data, created = PdbData.objects.get_or_create(pdb=rotamer_pdb)
                                    rotamer, created = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)

                                    fragment_data, created = PdbData.objects.get_or_create(pdb=fragment_pdb) 
                                    fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=fragment_data)
                                else:
                                    quit("Could not find "+f)

                                interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug=interactiontype,name=interactiontype)
                                fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure=structure,residue=residue,ligand=ligand,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)






            #print(results)
            results = sorted(results,key=itemgetter(3), reverse=True)



            #onlyfiles = [ f.replace('.txt','').split("_") for f in listdir(mypath) if isfile(join(mypath,f)) ]



            return render(request,'interaction/calculate.html',{'result' : "Looking at "+pdbname, 'outputs' : results })

        else:
            return HttpResponse("Error with form")
    else:
        return HttpResponse("Ooops how did you get here?")


def download(request):      
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')

    filename = module_dir+'/temp/results/'+pdbname+'/interaction/'+pdbname+'_'+ligand+'.pdb' # Select your file here.                                
    #wrapper = FileWrapper(file(filename))
    wrapper = FileWrapper(open(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Length'] = path.getsize(filename)
    return response

def pdb(request):       
    pdbname = request.GET.get('pdb')
    # response = HttpResponse(mimetype='application/force-download')
    # #response['Content-Disposition'] = 'attachment; filename=%s' % smart_str(file_name)
    # response['Content-Disposition'] = 'attachment; filename=%s' % smart_str(pdbname+'_'+ligand+'.pdb')
    # mypath = module_dir+'/temp/results/'+pdbname+'/interaction/'+pdbname+'_'+ligand+'.pdb'
    # response['X-Sendfile'] = smart_str(mypath)

    filename = module_dir+'/temp/pdbs/'+pdbname+'.pdb' # Select your file here.                                
    #wrapper = FileWrapper(file(filename))
    wrapper = FileWrapper(open(filename, 'rb'))
    response = HttpResponse(wrapper, content_type='text/plain')
    response['Content-Length'] = path.getsize(filename)
    return response

    # It's usually a good idea to set the 'Content-Length' header too.
    # You can also set any other required headers: Cache-Control, etc.
    #return response

    #context = {}
    #return HttpResponse(ref_id)
    #return render(request,'interaction/index.html',context)
