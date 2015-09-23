from django.shortcuts import render
from django.http import HttpResponse
#from interaction.functions import *
# Create your views here.
from subprocess import call
from .forms import PDBform
from django import forms
from django.core.servers.basehttp import FileWrapper
from django.db.models import Count, Min, Sum, Avg,Q


from interaction.models import *
from ligand.models import Ligand
from ligand.models import LigandType
from ligand.models import LigandRole, LigandProperities
from structure.models import Structure,PdbData,Rotamer,Fragment
from structure.functions import BlastSearch
from structure.assign_generic_numbers_gpcr import GenericNumbering
from protein.models import ProteinConformation,Protein
from residue.models import Residue, ResidueGenericNumber
from common.models import WebResource
from common.models import WebLink
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot

from os import path, listdir, devnull,makedirs
from os.path import isfile, join
import yaml
from operator import itemgetter
from datetime import datetime
import re
import json
import logging
from subprocess import Popen, DEVNULL
import urllib
import collections
from io import StringIO, BytesIO
from Bio.PDB import PDBIO, PDBParser

AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

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

def index(request, vignir=None):
    print(vignir)
    form = PDBform()

    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index','structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate( num_ligands=Count('structure_ligand_pair', distinct = True),num_interactions=Count('pk', distinct = True)).order_by('structure_ligand_pair__structure__pdb_code__index')
   
    #context = {}
    return render(request,'interaction/index.html',{'form': form, 'structures':structures, 'vignir':vignir })

def list(request):
    form = PDBform()
    #structures = ResidueFragmentInteraction.objects.distinct('structure_ligand_pair__structure').all()
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index','structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate( num_ligands=Count('structure_ligand_pair', distinct = True),num_interactions=Count('pk', distinct = True)).order_by('structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name')
    #structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index','structure_ligand_pair__structure').annotate(Count('structure_ligand_pair__ligand'))
    #print(structures.count())
    genes = {}
    countligands = {}
    totalligands = 0
    totalinteractions = 0
    totaltopinteractions = 0
    for structure in structures:
        #print(structure)
        if structure['structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name'] not in genes:
            genes[structure['structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name']] = 1
        totalligands += structure['num_ligands']
        totalinteractions += structure['num_interactions']
        ligands = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name').filter(structure_ligand_pair__structure__pdb_code__index=structure['structure_ligand_pair__structure__pdb_code__index']).annotate(numRes = Count('pk', distinct = True)).order_by('-numRes')
        for ligand in ligands:
            totaltopinteractions += ligand['numRes']
            if ligand['structure_ligand_pair__ligand__name'] not in countligands:
                countligands[ligand['structure_ligand_pair__ligand__name']] = 1
            break

        #print(structure.structure_ligand_pair.structure.pdb_code.index)
        #print(structure.numRes)
    #objects = Model.objects.filter(id__in=object_ids)
    #context = {}
    print('Structures with ligand information:' + str(structures.count()))
    print('Distinct genes:' + str(len(genes)))
    #print('ligands:' + str(totalligands))
    print('interactions:' + str(totalinteractions))
    print('interactions from top ligands:' + str(totaltopinteractions))
    print('Distinct ligands as top ligand:' + str(len(countligands)))

    return render(request,'interaction/list.html',{'form': form, 'structures': structures})

def crystal(request):
    pdbname = request.GET.get('pdb')
    form = PDBform()
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name').filter(structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes = Count('pk', distinct = True)).order_by('-numRes')
    crystal = Structure.objects.get(pdb_code__index=pdbname)
    p = Protein.objects.get(protein=crystal.protein_conformation.protein)
    residues = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname).order_by('rotamer__residue__sequence_number')
    print("residues",residues)
    return render(request,'interaction/crystal.html',{'form': form, 'pdbname': pdbname, 'structures': structures, 'crystal': crystal, 'protein':p, 'residues':residues })


def view(request):
    pdbname = request.GET.get('pdb')
    form = PDBform()
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name').filter(structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes = Count('pk', distinct = True)).order_by('-numRes')
    return render(request,'interaction/view.html',{'form': form, 'pdbname': pdbname, 'structures': structures})

def ligand(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    form = PDBform()
    fragments = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname).filter(structure_ligand_pair__ligand__name=ligand).order_by('interaction_type')
    return render(request,'interaction/ligand.html',{'form': form, 'pdbname': pdbname, 'ligand': ligand, 'fragments': fragments})

def fragment(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    fragment = request.GET.get('fragment')
    form = PDBform()
    fragments = ResidueFragmentInteraction.objects.get(id=fragment)
    return render(request,'interaction/fragment.html',{'form': form, 'pdbname': pdbname, 'ligand': ligand, 'fragmentid': fragment, 'fragments': fragments})

def updateall(request):
    structures = Structure.objects.values('pdb_code__index').distinct()
    for s in structures:
        pdbname = s['pdb_code__index']
        check = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname).all()

        if check.count()==0:
            t1 = datetime.now()
            runcalculation(pdbname)
            t2 = datetime.now()
            delta = t2 - t1
            seconds = delta.total_seconds()
            print("Calculation: Total time "+str(seconds)+" seconds for "+pdbname)
            t1 = datetime.now()
            results = parsecalculation(pdbname,False)
            t2 = datetime.now()
            delta = t2 - t1
            seconds = delta.total_seconds()
            print("Parsing: Total time "+str(seconds)+" seconds for "+pdbname)
            check = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname).all()
            print("Interactions found: " + str(check.count()))
        else:
            print(pdbname + " already calculated")
       
    #return render(request,'interaction/view.html',{'form': form, 'pdbname': pdbname, 'structures': structures})

def runcalculation(pdbname):
    call(["python", "interaction/functions.py","-p",pdbname], stdout=open(devnull, 'wb'), stderr=open(devnull, 'wb'))
    return None

def parsecalculation(pdbname, debug = True, ignore_ligand_preset = False): #consider skipping non hetsym ligands FIXME
    logger = logging.getLogger('build')
    mypath = '/tmp/interactions/results/'+pdbname+'/output'
    module_dir = '/tmp/interactions'
    results = []
    web_resource, created = WebResource.objects.get_or_create(slug='pdb',url='http://www.rcsb.org/pdb/explore/explore.do?structureId=$index')
    web_link, created = WebLink.objects.get_or_create(web_resource=web_resource,index=pdbname)

    structure=Structure.objects.filter(pdb_code=web_link) 
    if structure.exists():
        structure=Structure.objects.get(pdb_code=web_link)

        #quit() #quit!

        if structure.pdb_data is None:
            f = module_dir+"/pdbs/"+pdbname+".pdb"
            if isfile(f):      
                pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read()) #does this close the file?
            else:
                print('quitting due to no pdb in filesystem')
                quit()
            structure.pdb_data = pdbdata
            structure.save()

        protein=structure.protein_conformation

        for f in listdir(mypath):
            if isfile(join(mypath,f)):       
                result = yaml.load(open(mypath+"/"+f, 'rb'))
                output = result

                temp = f.replace('.yaml','').split("_")
                #print(output)
                temp.append([output])
                temp.append(round(output['score'][0][0]))
                temp.append((output['inchikey']).strip())
                temp.append((output['smiles']).strip())
                results.append(temp)

                if 'prettyname' not in output:
                    output['prettyname'] = temp[1]
                    #continue

                #print(' start ligand ' + output['prettyname'])
                ligand = Ligand.objects.filter(properities__inchikey=output['inchikey'].strip(), canonical=True)
                if ligand.exists():
                    ligand = ligand.get()
                    if output['prettyname']!=ligand.name: #add alias if same inchikey but different name.
                        alias = Ligand.objects.filter(name=output['prettyname'], properities__inchikey=output['inchikey'].strip(), canonical=False)
                        if alias.exists():
                            ligand = alias.get()
                        else:
                            alias = Ligand()
                            alias.name = output['prettyname']
                            alias.canonical = False
                            alias.properities = ligand.properities
                            alias.save()

                            ligand = alias #Use alias for structureligandinteraction
                else: #Ligand does not exist, create it
                    ligandtype, created = LigandType.objects.get_or_create(slug="sm", name='Small molecule') #FIXME
                    lp = LigandProperities()
                    lp.inchikey = output['inchikey'].strip()
                    lp.smiles = output['smiles'].strip()
                    lp.ligandtype = ligandtype
                    lp.save()

                    ligand = Ligand()
                    ligand.properities = lp
                    ligand.name = output['prettyname']
                    ligand.canonical = True #assume it's canonical but check.
                    ligand.ambigious_alias = False #assume till proven otherwise
                    ligand.save()
                    ligand.load_by_name(output['prettyname'])
                    ligand.save()

                    #if ligand.canonical== False: 
                        #print('looking for '+output['inchikey'].strip())
                        #ligand = Ligand.objects.get(properities__inchikey=output['inchikey'].strip(), canonical=True)
             
                #proteinligand, created = ProteinLigandInteraction.objects.get_or_create(protein=protein,ligand=ligand)

                f = module_dir+"/results/"+pdbname+"/interaction"+"/"+pdbname+"_"+temp[1]+".pdb"
                if isfile(f):      
                    pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read()) #does this close the file?
                    if debug: print("Found file"+f)
                else:
                    print('quitting due to no pdb for fragment in filesystem',f)
                    quit()


                structureligandinteraction = StructureLigandInteraction.objects.filter(ligand__properities__inchikey=output['inchikey'].strip(),structure=structure)
                if structureligandinteraction.exists():
                    structureligandinteraction = structureligandinteraction.get()
                    structureligandinteraction.pdb_file = pdbdata
                    structureligandinteraction.pdb_reference = temp[1]
                elif StructureLigandInteraction.objects.filter(pdb_reference=temp[1],structure=structure).exists(): #incase defined reference doesn't match on inchikey
                    structureligandinteraction = StructureLigandInteraction.objects.filter(pdb_reference=temp[1],structure=structure).get()
                    structureligandinteraction.pdb_file = pdbdata
                    if structureligandinteraction.ligand.properities.inchikey is None:
                        logger.info('Old ligand didnt get inchikey -- error in naming, using inchikey/properities from structure')
                        structureligandinteraction.ligand.delete()
                        structureligandinteraction.ligand = ligand
                else:
                    ligandrole, created = LigandRole.objects.get_or_create(name='unknown',slug='unknown')
                    structureligandinteraction = StructureLigandInteraction()
                    structureligandinteraction.ligand = ligand
                    structureligandinteraction.structure = structure
                    structureligandinteraction.ligand_role = ligandrole
                    structureligandinteraction.pdb_file = pdbdata
                    structureligandinteraction.pdb_reference = temp[1]

                structureligandinteraction.save()
                

                #structureligandinteraction, created = StructureLigandInteraction.objects.get_or_create(ligand=ligand,structure=structure, ligand_role=ligandrole, pdb_file=pdbdata) #, pdb_reference=pdbname <-- max length set to 3? So can't insert ones correctly

                
                for interactiontype,interactionlist in output.items():
                    if interactiontype=='hbond' or interactiontype=='hbondplus':
                        for entry in interactionlist:
                            #print(interactiontype,entry)
                            aa = entry[0]
                            aa,pos,chain = regexaa(aa)

                            residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos)
                            if residue.exists():
                                residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos)
                                if residue.amino_acid!=aa:
                                    if debug: print("Updated amino acid from",residue.amino_acid,"to",aa)
                                    residue.amino_acid = aa
                                    residue.save()
                            else:
                                if debug: print("Not found residue!",pdbname,pos,aa)
                                residue, created=Residue.objects.get_or_create(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                #continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                            for pair in entry[3]:
                                fragment = pair[0] #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                                f = module_dir+"/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_"+fragment+"_HB.pdb"
                                if interactiontype=='hbondplus':  f = module_dir+"/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_"+fragment+"_HBC.pdb"
                                if isfile(f):      
                                    if debug: print("Found file"+f)
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
                                    fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=fragment_data, residue=residue)
                                else:
                                    quit("Could not find "+f)

                                interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug=interactiontype,name=interactiontype)
                                fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure_ligand_pair=structureligandinteraction,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)
                                
                                
                    elif interactiontype=='hbond_confirmed':
                        for entry in interactionlist:
                            #print(interactiontype,entry)
                            aa = entry[0]
                            aa,pos,chain = regexaa(aa)
                            interactiontype="HB"+entry[1][0][0]

                            residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos)
                            if residue.exists():
                                residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos)
                                if residue.amino_acid!=aa:
                                    if debug: logger.info("Updated amino acid from",residue.amino_acid,"to",aa)
                                    residue.amino_acid = aa
                                    residue.save()
                            else:
                                if debug: logger.info("Not found residue!",pdbname,pos,aa)
                                residue, created=Residue.objects.get_or_create(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                #continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                            for pair in entry[1]:  
                                fragment = pair[1] #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                                #mypath = module_dir+'/temp/results/'+pdbname+'/fragments'
                                f = module_dir+"/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_"+fragment+"_HB.pdb"
                                if isfile(f):      
                                    if debug: print("Found file"+f)
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
                                    fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=fragment_data, residue=residue)
                                else:
                                    quit("Could not find "+f)

                                interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug=interactiontype,name=interactiontype)
                                #fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure=structure,residue=residue,ligand=ligand,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)
                                fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure_ligand_pair=structureligandinteraction,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)
                    elif interactiontype=='hydrophobic': #NO FRAGMENT FOR THESE, WHAT TO DO? USE WHOLE LIGAND?
                        for entry in interactionlist:
                            if debug: print(entry)
                            aa = entry[0]
                            aa,pos,chain = regexaa(aa)
                            #fragment = entry[1][0][1] #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                            residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos)
                            if residue.exists():
                                residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos)
                                if residue.amino_acid!=aa:
                                    if debug: logger.info("Updated amino acid from",residue.amino_acid,"to",aa)
                                    residue.amino_acid = aa
                                    residue.save()
                            else:
                                if debug: logger.info("Not found residue!",pdbname,pos,aa)
                                residue, created=Residue.objects.get_or_create(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                #continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                            fragment = '' #NEED TO EXPAND THIS TO INCLUDE MORE INFO

                            f = module_dir+"/results/"+pdbname+"/ligand/"+temp[1]+"_"+pdbname+".pdb"
                            if isfile(f):      
                                liganddata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read()) #does this close the file?
                                if debug: logger.info("Hydro Found file"+f)
                            else:
                                quit()

                            f = module_dir+"/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"__hydrop.pdb"
                            rotamer_data, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read())

                            rotamer, created = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)

                            fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=liganddata, residue=residue)

                            interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug='hydrofob',name=interactiontype)
                            #fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure=structure,residue=residue,ligand=ligand,interaction_type=interaction_type,fragment=fragment)
                            fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure_ligand_pair=structureligandinteraction,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)
                    elif interactiontype=='aromaticplus' or interactiontype=='aromatic' or interactiontype=='aromaticfe':
                        for entry in interactionlist:
                            if debug: logger.info(entry)
                            aa = entry[0]
                            aa,pos,chain = regexaa(aa)
                            fragment = entry[1]

                            residue=Residue.objects.filter(protein_conformation=protein, sequence_number=pos)
                            if residue.exists():
                                residue=Residue.objects.get(protein_conformation=protein, sequence_number=pos)
                                if residue.amino_acid!=aa:
                                    if debug: logger.info("Updated amino acid from",residue.amino_acid,"to",aa)
                                    residue.amino_acid = aa
                                    residue.save()
                            else:
                                if debug: logger.info("Not found residue!",pdbname,pos,aa)
                                residue, created=Residue.objects.get_or_create(protein_conformation=protein, sequence_number=pos,amino_acid=aa)
                                #continue #SKIP THESE -- mostly fusion residues that aren't mapped yet.

                            f = module_dir+"/results/"+pdbname+"/fragments"+"/"+pdbname+"_"+temp[1]+"_"+entry[0]+"_aromatic_"+str(entry[1])+".pdb"
                            if isfile(f):      
                                if debug: logger.info("Found file"+f)
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
                                fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=fragment_data, residue=residue)
                            else:
                                quit("Could not find "+f)

                            interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(slug=interactiontype,name=interactiontype)
                            #fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure=structure,residue=residue,ligand=ligand,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)
                            fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(structure_ligand_pair=structureligandinteraction,interaction_type=interaction_type,fragment=fragment, rotamer=rotamer)

    else:
        if debug: logger.info("Structure not in DB?!??!")
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



            #print(results)
    results = sorted(results,key=itemgetter(3), reverse=True)

    return results

def runusercalculation(filename,session):
    call(["python", "interaction/functions.py","-p",filename,"-s",session])
    return None

def parseusercalculation(pdbname,session, debug = True, ignore_ligand_preset = False, ): #consider skipping non hetsym ligands FIXME
    logger = logging.getLogger('build')
    mypath = '/tmp/interactions/'+session+'/results/'+pdbname+'/output'
    module_dir = '/tmp/interactions/'+session
    results = []
   

    for f in listdir(mypath):
        if isfile(join(mypath,f)):       
            result = yaml.load(open(mypath+"/"+f, 'rb'))
            output = result

            temp = f.replace('.yaml','').split("_")
            #print(output)
            temp.append([output])
            temp.append(round(output['score'][0][0]))
            temp.append((output['inchikey']).strip())
            temp.append((output['smiles']).strip())
            results.append(temp)

            if 'prettyname' not in output:
                output['prettyname'] = temp[1]
                #continue

            #print(' start ligand ' + output['prettyname'])
        
            #print(results)
    results = sorted(results,key=itemgetter(3), reverse=True)

    return results

def calculate(request, vignir=None):   
    if request.method == 'POST':
        form = PDBform(request.POST, request.FILES)
        if form.is_valid():

            pdbname = form.cleaned_data['pdbname'].strip()
            results = ''

            module_dir = '/tmp/interactions/' + request.session.session_key

            if not path.exists('/tmp/interactions/'):
                makedirs('/tmp/interactions/')
            if not path.exists('/tmp/interactions/'+ request.session.session_key):
                makedirs('/tmp/interactions/'+ request.session.session_key)
            if not path.exists('/tmp/interactions/'+ request.session.session_key+'/pdbs'):
                makedirs('/tmp/interactions/'+ request.session.session_key+'/pdbs')
            if not path.exists('/tmp/interactions/'+ request.session.session_key+'/temp'):
                makedirs('/tmp/interactions/'+ request.session.session_key+'/temp')

            if 'file' in request.FILES:
                pdbdata = request.FILES['file']
                pdbname = path.splitext(str(pdbdata))[0]
                with open(module_dir+'/pdbs/'+str(pdbdata), 'wb+') as destination:
                     for chunk in pdbdata.chunks():
                         destination.write(chunk)

                print(module_dir+'/pdbs/'+str(pdbdata))
                temp_path = module_dir+'/pdbs/'+str(pdbdata)
                pdbdata=open(temp_path,'r').read()
                print('do calc')
                runusercalculation(pdbname,request.session.session_key)

            else:
                pdbname = form.cleaned_data['pdbname'].strip()
                print('pdbname selected!',pdbname)
                temp_path = module_dir+'/pdbs/'+pdbname+'.pdb'

                if not path.isfile(temp_path):
                    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbname
                    pdbdata = urllib.request.urlopen(url).read().decode('utf-8')
                    f=open(temp_path,'w')
                    f.write(pdbdata)
                    f.close();
                else:
                    pdbdata=open(temp_path,'r').read()
                print('do calc')
                runusercalculation(pdbname,request.session.session_key)

            print('map GPCRdb numbering on result')
            #MAPPING GPCRdb numbering onto pdb.
            generic_numbering = GenericNumbering(temp_path)
            out_struct = generic_numbering.assign_generic_numbers()
            structure_residues = generic_numbering.residues
            segments = {}

            generic_ids = []
            generic_number = []
            previous_seg = 'N-term'

            # for chain in out_struct:
            #     for residue in chain:
            #         if residue.id[1] in structure_residues[chain.id].keys():
            #             print(residue)
            #         else:
            #             print('not found')
            #             print(residue)

            for c,res in structure_residues.items():
                for i,r in sorted(res.items()): #sort to be able to assign loops
                    if r.gpcrdb:
                        if r.gpcrdb[0]=='-':
                            r.gpcrdb = r.gpcrdb[1:]+"1" #fix stefan - for bulge
                        r.gpcrdb = str(r.gpcrdb).replace('.','x')
                        generic_number.append(r.gpcrdb)
                    if r.gpcrdb_id:
                        generic_ids.append(r.gpcrdb_id)
                    if (r.segment):
                        if not r.segment in segments: 
                            segments[r.segment] = {}
                        segments[r.segment][r.number] = [r.display,r.name,r.gpcrdb]
                        previous_seg = r.segment
                    else: #if no segment assigned by blast
                        if previous_seg in ['N-term','ICL1','ECL1','ICL2','ECL2','ICL3','ICL3','C-term']:
                            if not previous_seg in segments: 
                                segments[previous_seg] = {}
                            segments[previous_seg][r.number] = ['',r.name, '']
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
                            segments[previous_seg][r.number] = ['',r.name, '']

            #print(segments)
            #print(generic_ids)
            rs = ResidueGenericNumber.objects.filter(id__in=generic_ids)

            reference_generic_numbers = {}
            for r in rs: #make lookup dic.
                reference_generic_numbers[r.label] = r

            rs = ResidueGenericNumber.objects.filter(label__in=generic_number)
            reference_numbers = {}
            for r in rs: #make lookup dic.
                reference_numbers[r.label] = r
            ################################################################################
            #print(reference_numbers)
            residue_list=[]
            for seg,reslist in segments.items():

                for seq_number,v in sorted(reslist.items()):
                    r = Residue()
                    r.sequence_number =  seq_number
                    if "x" in v[0]:
                        r.display_number = reference_numbers[v[2]] #FIXME
                        r.display_generic_number = reference_generic_numbers[v[0]] #FIXME
                        r.segment_slug = seg
                        r.family_generic_number = v[2]
                    else:
                        r.segment_slug = seg
                    r.amino_acid = v[1]
                    residue_list.append(r)
     
            #print(residue_list)

            HelixBox = DrawHelixBox(residue_list,'Class A',str('test'), nobuttons = 1)
            SnakePlot = DrawSnakePlot(residue_list,'Class A',str('test'), nobuttons = 1)

            xtal = {}
            hetsyn = {}
            hetsyn_reverse = {}
            for line in pdbdata.splitlines():
                if line.startswith('HETSYN'): 
                    m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)",line) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
                    if (m):
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('HETNAM'): 
                    m = re.match("HETNAM[\s]+([\w]{3})[\s]+(.+)",line) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
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

            #print(xtal)

            results = parseusercalculation(pdbname,request.session.session_key)

            #print(results)
            simple = collections.OrderedDict()
            simple_generic_number = collections.OrderedDict()
            residues = []
            mainligand = ''
            for ligand in results:
                print(ligand[1])
                if mainligand=='': mainligand=ligand[1] #select top hit
                simple[ligand[1]] = {'score':round(ligand[2][0]['score'][0][0])}
                simple_generic_number[ligand[1]] = {'score':round(ligand[2][0]['score'][0][0])}
                for key,values in ligand[2][0].items():
                    if key in ['aromatic','aromaticplus','hbond','hbond_confirmed','hydrophobic', 'hbondplus', 'aromaticfe','waals']:
                        for value in values:
                            aa,pos,chain = regexaa(value[0])
                            if int(pos) in structure_residues[chain]:
                                r = structure_residues[chain][int(pos)]
                                display = r.display
                                segment = r.segment

                                if display!="" and display in simple_generic_number[ligand[1]]:
                                    simple_generic_number[ligand[1]][display].append(key)
                                elif display!="":
                                    simple_generic_number[ligand[1]][display] = [key]
                            else:
                                display = ''
                                segment = ''

                            if value[0] in simple[ligand[1]]:
                                simple[ligand[1]][value[0]].append(key)
                            else:
                                simple[ligand[1]][value[0]] = [key]

                            residues.append({'type':key,'aa':aa,'ligand':ligand[1],'pos':pos, 'gpcrdb':display, 'segment':segment})

            if vignir:
                return render(request,'interaction/vignir.html',{'simple':simple,'simple_generic_number':simple_generic_number})
            else:
                return render(request,'interaction/diagram.html',{'result' : "Looking at "+pdbname, 'outputs' : results,
                 'simple' : simple ,'simple_generic_number' : simple_generic_number , 'xtal' : xtal, 'pdbname':pdbname, 'mainligand':mainligand, 'residues':residues,
                 'HelixBox':HelixBox, 'SnakePlot':SnakePlot})

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
        pdbdata = open('/tmp/interactions/'+session+'/results/'+pdbname+'/interaction/'+pdbname+'_'+ligand+'.pdb', 'r').read()
        response=HttpResponse(pdbdata, content_type='text/plain')
    else:

        pair = StructureLigandInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(Q(ligand__properities__inchikey=ligand) | Q(ligand__name=ligand)).exclude(pdb_file__isnull=True).get()
        response = HttpResponse(pair.pdb_file.pdb, content_type='text/plain')
    return response

def ajax(request, slug, **response_kwargs):
    interactions = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name=slug).order_by('rotamer__residue__sequence_number')
    print(interactions)
    #return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for interaction in interactions:
        sequence_number = interaction.rotamer.residue.sequence_number
        aa = interaction.rotamer.residue.amino_acid
        interactiontype = interaction.interaction_type.name
        if sequence_number not in jsondata: jsondata[sequence_number] = []
        jsondata[sequence_number].append([aa,interactiontype])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def ajaxLigand(request, slug, ligand, **response_kwargs):
    print(ligand)
    interactions = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name=slug,structure_ligand_pair__ligand__name=ligand).order_by('rotamer__residue__sequence_number')
    print(interactions)
    #return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for interaction in interactions:
        sequence_number = interaction.rotamer.residue.sequence_number
        aa = interaction.rotamer.residue.amino_acid
        interactiontype = interaction.interaction_type.name
        if sequence_number not in jsondata: jsondata[sequence_number] = []
        jsondata[sequence_number].append([aa,interactiontype])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def pdbfragment(request):      
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    fragment = request.GET.get('fragment')

    result = ResidueFragmentInteraction.objects.filter(id=fragment).get()
    response = HttpResponse(result.rotamer.pdbdata.pdb+result.fragment.pdbdata.pdb, content_type='text/plain')
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
        pdbdata = open('/tmp/interactions/'+session+'/pdbs/'+pdbname+'.pdb', 'r').read()
        response=HttpResponse(pdbdata, content_type='text/plain')
    else:
        web_resource, created = WebResource.objects.get_or_create(slug='pdb',url='http://www.rcsb.org/pdb/explore/explore.do?structureId=$index')
        web_link, created = WebLink.objects.get_or_create(web_resource=web_resource,index=pdbname)

        structure=Structure.objects.filter(pdb_code=web_link) 
        if structure.exists():
            structure=Structure.objects.get(pdb_code=web_link)
        else:
             quit() #quit!

        if structure.pdb_data is None:
            quit()

        response = HttpResponse(structure.pdb_data.pdb, content_type='text/plain')
    return response
