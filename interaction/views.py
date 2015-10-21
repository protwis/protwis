from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django import forms
from django.core.servers.basehttp import FileWrapper
from django.db.models import Count, Min, Sum, Avg,Q
from django.utils.text import slugify


from interaction.models import *
from interaction.forms import PDBform
from ligand.models import Ligand
from ligand.models import LigandType
from ligand.models import LigandRole, LigandProperities
from structure.models import Structure,PdbData,Rotamer,Fragment
from structure.functions import BlastSearch
from structure.assign_generic_numbers_gpcr import GenericNumbering
from protein.models import ProteinConformation,Protein, ProteinSegment
from residue.models import Residue, ResidueGenericNumber, ResidueGenericNumberEquivalent,ResidueNumberingScheme
from common.models import WebResource
from common.models import WebLink
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot
from common.selection import SimpleSelection, Selection, SelectionItem
from common import definitions

import os
from os import listdir, devnull, makedirs
import yaml
from operator import itemgetter
from datetime import datetime
import re
import json
import logging
from subprocess import call, Popen, DEVNULL
import urllib
import collections
from collections import OrderedDict
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

def index(request):
    form = PDBform()

    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index','structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate( num_ligands=Count('structure_ligand_pair', distinct = True),num_interactions=Count('pk', distinct = True)).order_by('structure_ligand_pair__structure__pdb_code__index')
   
    #context = {}
    return render(request,'interaction/index.html', {'form': form, 'structures':structures})

def list_structures(request):
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
    calc_script = os.sep.join([os.path.dirname(__file__), 'functions.py'])
    call(["python", calc_script, "-p",pdbname], stdout=open(devnull, 'wb'), stderr=open(devnull, 'wb'))
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
            if os.path.isfile(f):      
                pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read()) #does this close the file?
            else:
                print('quitting due to no pdb in filesystem')
                quit()
            structure.pdb_data = pdbdata
            structure.save()

        protein=structure.protein_conformation

        for f in listdir(mypath):
            if os.path.isfile(os.path.join(mypath,f)):       
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
                    output['prettyname'] = temp[1] #use hetsyn name if possible, others 3letter
                    #continue

                f = module_dir+"/results/"+pdbname+"/interaction"+"/"+pdbname+"_"+temp[1]+".pdb"
                if os.path.isfile(f):      
                    pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read()) #does this close the file?
                    if debug: print("Found file"+f)
                else:
                    print('quitting due to no pdb for fragment in filesystem',f)
                    quit()

                structureligandinteraction = StructureLigandInteraction.objects.filter(pdb_reference=temp[1],structure=structure, annotated=True)
                if structureligandinteraction.exists(): #if the annotated exists
                    structureligandinteraction = structureligandinteraction.get()
                    structureligandinteraction.pdb_file = pdbdata
                    ligand = structureligandinteraction.ligand
                    if structureligandinteraction.ligand.properities.inchikey!=output['inchikey'].strip():
                        logger.error('inchikey for annotated ligand and PDB ligand mismatch' + output['prettyname'])
                elif StructureLigandInteraction.objects.filter(pdb_reference=temp[1],structure=structure).exists():
                    structureligandinteraction = StructureLigandInteraction.objects.filter(pdb_reference=temp[1],structure=structure).get()
                    structureligandinteraction.pdb_file = pdbdata
                    ligand = structureligandinteraction.ligand
                else: #create ligand and pair

                    ligand = Ligand.objects.filter(name=output['prettyname'], canonical=True)

                    if ligand.exists(): #if ligand with name (either hetsyn or 3 letter) exists use that.
                        ligand = ligand.get()
                    else: #create it
                        default_ligand_type = 'N/A'
                        lt, created = LigandType.objects.get_or_create(slug=slugify(default_ligand_type),
                            defaults={'name': default_ligand_type})

                        ligand = Ligand()
                        ligand = ligand.load_from_pubchem('inchikey', output['inchikey'].strip(), lt, output['prettyname'])
                        ligand.save()

                    ligandrole, created = LigandRole.objects.get_or_create(name='unknown',slug='unknown')
                    structureligandinteraction = StructureLigandInteraction()
                    structureligandinteraction.ligand = ligand
                    structureligandinteraction.structure = structure
                    structureligandinteraction.ligand_role = ligandrole
                    structureligandinteraction.pdb_file = pdbdata
                    structureligandinteraction.pdb_reference = temp[1]

                structureligandinteraction.save()
                
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
                                if os.path.isfile(f):      
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
                                if os.path.isfile(f):      
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
                            if os.path.isfile(f):      
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
                            if os.path.isfile(f):      
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
            if os.path.isfile(os.path.join(mypath,f)):       
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

def runusercalculation(filename, session):
    calc_script = os.sep.join([os.path.dirname(__file__), 'functions.py'])
    call(["python", calc_script,"-p",filename,"-s",session])
    return None

def parseusercalculation(pdbname, session, debug = True, ignore_ligand_preset = False, ): #consider skipping non hetsym ligands FIXME
    logger = logging.getLogger('build')
    mypath = '/tmp/interactions/'+session+'/results/'+pdbname+'/output'
    module_dir = '/tmp/interactions/'+session
    results = []
   

    for f in listdir(mypath):
        if os.path.isfile(os.path.join(mypath,f)):       
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

def calculate(request, redirect=None):   
    if request.method == 'POST':
        form = PDBform(request.POST, request.FILES)
        if form.is_valid():

            pdbname = form.cleaned_data['pdbname'].strip()
            results = ''

            if not request.session.exists(request.session.session_key):
                request.session.create() 
            session_key = request.session.session_key
            
            module_dirs = []
            module_dir = '/tmp/interactions'
            module_dirs.append(module_dir)
            module_dir = os.sep.join([module_dir, session_key])
            module_dirs.append(module_dir)
            module_dirs.append(os.sep.join([module_dir, 'pdbs']))
            module_dirs.append(os.sep.join([module_dir, 'temp']))
            
            # create dirs and set permissions (needed on some systems)
            for mdir in module_dirs:
                os.makedirs(mdir, exist_ok=True)
                os.chmod(mdir, 0o777)

            if 'file' in request.FILES:
                pdbdata = request.FILES['file']
                pdbname = os.path.splitext(str(pdbdata))[0]
                with open(module_dir+'/pdbs/'+str(pdbdata), 'wb+') as destination:
                     for chunk in pdbdata.chunks():
                         destination.write(chunk)

                temp_path = module_dir+'/pdbs/'+str(pdbdata)
                pdbdata=open(temp_path,'r').read()
                runusercalculation(pdbname, session_key)

            else:
                pdbname = form.cleaned_data['pdbname'].strip()
                #print('pdbname selected!',pdbname)
                temp_path = module_dir+'/pdbs/'+pdbname+'.pdb'

                if not os.path.isfile(temp_path):
                    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbname
                    pdbdata = urllib.request.urlopen(url).read().decode('utf-8')
                    f=open(temp_path,'w')
                    f.write(pdbdata)
                    f.close();
                else:
                    pdbdata=open(temp_path,'r').read()
                runusercalculation(pdbname, session_key)

            # MAPPING GPCRdb numbering onto pdb.
            generic_numbering = GenericNumbering(temp_path)
            out_struct = generic_numbering.assign_generic_numbers()
            structure_residues = generic_numbering.residues
            prot_id_list = generic_numbering.prot_id_list
            segments = {}

            generic_ids = []
            generic_number = []
            previous_seg = 'N-term'

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

            rs = ResidueGenericNumber.objects.filter(id__in=generic_ids)

            reference_generic_numbers = {}
            for r in rs: #make lookup dic.
                reference_generic_numbers[r.label] = r

            rs = ResidueGenericNumber.objects.filter(label__in=generic_number)
            reference_numbers = {}
            for r in rs: #make lookup dic.
                reference_numbers[r.label] = r

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
                if line.startswith('REMARK   2 RESOLUTION.'):
                    xtal['resolution'] = line[22:].strip()

            results = parseusercalculation(pdbname,request.session.session_key)

            simple = collections.OrderedDict()
            simple_generic_number = collections.OrderedDict()
            residues_browser = []
            residue_table_list = []
            mainligand = ''
            interaction_types = ['aromatic','aromaticplus','hbond','hbond_confirmed','hydrophobic', 'hbondplus',
                'aromaticfe','waals']
            
            for ligand in results:
                ligand_score = round(ligand[2][0]['score'][0][0])

                # select top hit
                if mainligand == '':
                    mainligand=ligand[1]
                
                simple[ligand[1]] = {'score': ligand_score}
                simple_generic_number[ligand[1]] = {'score': ligand_score}
                
                for key,values in ligand[2][0].items():
                    if key in interaction_types:
                        for value in values:
                            aa,pos,chain = regexaa(value[0])
                            if int(pos) in structure_residues[chain]:
                                r = structure_residues[chain][int(pos)]
                                display = r.display
                                segment = r.segment
                                generic = r.gpcrdb

                                if generic!="":
                                    residue_table_list.append(generic)

                                if generic!="" and generic in simple_generic_number[ligand[1]]:
                                    simple_generic_number[ligand[1]][generic].append(key)
                                elif generic!="":
                                    simple_generic_number[ligand[1]][generic] = [key]
                            else:
                                display = ''
                                segment = ''

                            if value[0] in simple[ligand[1]]:
                                simple[ligand[1]][value[0]].append(key)
                            else:
                                simple[ligand[1]][value[0]] = [key]

                            residues_browser.append({'type':key,'aa':aa,'ligand':ligand[1],'pos':pos, 'gpcrdb':display, 'segment':segment})
                break #only use the top one

            #RESIDUE TABLE
            segments = ProteinSegment.objects.all().filter().prefetch_related()
            #s = Structure.objects.get(pdb_code__index=xtal['pdb_code'])
            #proteins = [s.protein_conformation.protein]
            proteins = []
            protein_list = Protein.objects.filter(pk__in=prot_id_list)
            numbering_schemes_selection = [settings.DEFAULT_NUMBERING_SCHEME]
            for p in protein_list:
                proteins.append(p)
                if p.residue_numbering_scheme.slug not in numbering_schemes_selection:
                    numbering_schemes_selection.append(p.residue_numbering_scheme.slug)

            numbering_schemes = ResidueNumberingScheme.objects.filter(slug__in=numbering_schemes_selection).all()
            default_scheme = numbering_schemes.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
            data = OrderedDict()

            for segment in segments:
                data[segment.slug] = OrderedDict()
                residues = Residue.objects.filter(protein_segment=segment,  protein_conformation__protein__in=proteins, 
                    generic_number__label__in=residue_table_list).prefetch_related('protein_conformation__protein', 
                    'protein_conformation__state', 'protein_segment',
                    'generic_number','display_generic_number','generic_number__scheme',
                    'alternative_generic_numbers__scheme')
                for scheme in numbering_schemes:
                    if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                        for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                            data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}
                    elif scheme == default_scheme:
                        for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                            data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}

                for residue in residues:
                    alternatives = residue.alternative_generic_numbers.all()
                    pos = residue.generic_number
                    for alternative in alternatives:
                        if alternative.scheme not in numbering_schemes:
                            continue
                        scheme = alternative.scheme
                        if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                            pos = residue.generic_number
                            if scheme == pos.scheme:
                                data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                            else:
                                if scheme.slug not in data[segment.slug][pos.label].keys():
                                    data[segment.slug][pos.label][scheme.slug] = alternative.label
                                if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                    data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                                data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                        else:
                            if scheme.slug not in data[segment.slug][pos.label].keys():
                                data[segment.slug][pos.label][scheme.slug] = alternative.label
                            if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                            data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)

            # Preparing the dictionary of list of lists. Dealing with tripple nested dictionary in django templates is a nightmare
            flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])
            for s in iter(flattened_data):
                flattened_data[s] = [[data[s][x][y.slug] for y in numbering_schemes]+data[s][x]['seq'] for x in sorted(data[s])]
            
            context = {}
            context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [x.name for x in proteins],[x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
            context['segments'] = [x.slug for x in segments if len(data[x.slug])]
            context['data'] = flattened_data
            context['number_of_schemes'] = len(numbering_schemes)
   
            if redirect:
                # get simple selection from session
                simple_selection = request.session.get('selection', False)
                
                # create full selection and import simple selection (if it exists)
                selection = Selection()
                if simple_selection:
                    selection.importer(simple_selection)

                # convert identified interactions to the correct names and add them to the session
                # FIXME use the correct names from the beginning
                interaction_name_dict = {
                    'aromatic': 'ar',
                    'aromaticplus': 'ar',
                    'aromaticfe': 'ar',
                    'hbond': 'hba',
                    'hbond_confirmed': 'hba',
                    'hbondplus': 'hba',
                    # 'hydrophobic': 'hp',
                }
                for gn, interactions in simple_generic_number[mainligand].items():
                    print(gn, interactions)
                    if gn != 'score' and gn != 0.0: # FIXME leave these out when dict is created
                        for interaction in interactions:
                            if interaction in interaction_name_dict:
                                feature = interaction_name_dict[interaction]
                                break
                        else:
                            continue

                        # get residue number equivalent object
                        rne = ResidueGenericNumberEquivalent.objects.get(label=gn,
                            scheme__slug='gpcrdba')

                        # create a selection item
                        properties = {
                            'feature': feature,
                            'amino_acids': ','.join(definitions.AMINO_ACID_GROUPS[feature])
                        }
                        selection_item = SelectionItem('site_residue', rne, properties)

                        # add to selection
                        selection.add('segments', 'site_residue', selection_item)

                # export simple selection that can be serialized
                simple_selection = selection.exporter()

                # add simple selection to session
                request.session['selection'] = simple_selection
                
                # re-direct to segment selection (with the extracted interactions already selected)
                return HttpResponseRedirect(redirect)
            else:
                return render(request,'interaction/diagram.html',{'result' : "Looking at "+pdbname, 'outputs' : results,
                 'simple' : simple ,'simple_generic_number' : simple_generic_number , 'xtal' : xtal, 'pdbname':pdbname, 'mainligand':mainligand, 'residues':residues_browser,
                 'HelixBox':HelixBox, 'SnakePlot':SnakePlot ,'data':context['data'], 
                'header':context['header'], 'segments':context['segments'], 'number_of_schemes':len(numbering_schemes)})

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
