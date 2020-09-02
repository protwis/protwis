from django.utils.text import slugify
from django.db import IntegrityError

#from chembl_webresource_client import new_client
from common.models import WebResource
from common.models import WebLink
from ligand.models import Ligand, LigandType, LigandProperities

def get_or_make_ligand(ligand_id,type_id, name = None):
    if type_id=='PubChem CID' or type_id=='SMILES':
        if type_id=='PubChem CID':
            pubchem_lookup_value = 'cid'
        elif type_id=='SMILES':
            pubchem_lookup_value = 'smiles'

        try:
            web_resource = WebResource.objects.get(slug='pubchem')
        except:
            # abort if pdb resource is not found
            raise Exception('PubChem resource not found, aborting!')
        if name:
            ligand_name = name
        else:
            ligand_name = False

        try:
            # if this name is canonical and it has a ligand record already
            if (ligand_name==False):

                l = None
                ls = Ligand.objects.filter(canonical=True,
                   properities__web_links__web_resource=web_resource,
                   properities__web_links__index=ligand_id)

                for ligand in ls:
                    l = ligand
                    #print (l)
                    break
                if l == None:
                    l = Ligand.objects.get(canonical=True,
                    properities__web_links__web_resource=web_resource,
                    properities__web_links__index=ligand_id)

            else:
               l = Ligand.objects.get(name=ligand_name, canonical=True,
                   properities__web_links__web_resource=web_resource,
                   properities__web_links__index=ligand_id)

            #l = Ligand.objects.get(name=ligand_name, canonical=True,
            #    properities__web_links__web_resource=web_resource,
            #    properities__web_links__index=ligand_id)
            #
        except Ligand.DoesNotExist:
            try:
                # if exists under different name
                l_canonical = Ligand.objects.get(properities__web_links__web_resource=web_resource,
                    properities__web_links__index=ligand_id, canonical=True)
                #print (created)
                try:
                    l, created = Ligand.objects.get_or_create(properities = l_canonical.properities,
                        name = ligand_name, canonical = False)
                except IntegrityError:
                    l = Ligand.objects.get(properities = l_canonical.properities,
                        name = ligand_name, canonical = False)
            except Ligand.DoesNotExist:
                # fetch ligand from pubchem
                default_ligand_type = 'Small molecule'
                lt, created = LigandType.objects.get_or_create(slug=slugify(default_ligand_type),
                    defaults={'name': default_ligand_type})
                l = Ligand()
                #print (ligand_name)
                l = l.load_from_pubchem(pubchem_lookup_value, ligand_id, lt, ligand_name)
                #print (l)
                if l == None and type_id=='SMILES': #insert manually if smiles and unfound in pubchem
                    try:
                        l = Ligand.objects.get(name=ligand_name, canonical=True,
                                                properities__smiles=ligand_id)
                    except Ligand.DoesNotExist:
                        try:
                            l = Ligand.objects.get(name__startswith=ligand_name, canonical=True,properities__smiles=ligand_id) #if no properities exist
                        except Ligand.DoesNotExist:
                            try:
                                l = Ligand.objects.get(name=ligand_name, canonical=True,properities__smiles=None) #if no properities exist
                                l.properities.smiles = ligand_id
                                l.properities.save()
                                l.save()
                            except Ligand.DoesNotExist:
                                ## now insert a new ligand, but first make sure name is unique
                                if Ligand.objects.filter(name=ligand_name).exists():
                                    ls = Ligand.objects.filter(name__startswith=ligand_name, canonical=True).order_by("pk")
                                    for l_temp in ls:
                                        last = l_temp.name.split("_")[-1]
                                    if last==ligand_name: #no addition yet
                                        ligand_name = ligand_name +"_1"
                                    else:
                                        ligand_name = ligand_name +"_"+str(int(last)+1)
                                l = Ligand()
                                l.name = ligand_name
                                lp = LigandProperities()
                                lp.smiles = ligand_id
                                lp.ligand_type = lt
                                lp.save()
                                l.properities = lp
                                l.canonical = True #maybe false, but that would break stuff.
                                l.ambigious_alias = False
                                try:
                                    l.save()
                                except IntegrityError:
                                    l = Ligand.objects.get(name=ligand_name, canonical=True)

    elif name:

        # if this name is canonical and it has a ligand record already
        if Ligand.objects.filter(name=name, canonical=True).exists():
            l = Ligand.objects.get(name=name, canonical=True)

        # if this matches an alias that only has "one" parent canonical name - eg distinct
        elif Ligand.objects.filter(name=name, canonical=False,
            ambigious_alias=False).exists():
            l = Ligand.objects.get(name=name, canonical=False, ambigious_alias=False)

        # if this matches an alias that only has several canonical parents, must investigate, start
        # with empty.
        elif Ligand.objects.filter(name=name, canonical=False,
            ambigious_alias=True).exists():
            lp = LigandProperities()
            lp.save()
            l = Ligand()
            l.properities = lp
            l.name = name
            l.canonical = False
            l.ambigious_alias = True
            l.save()
            l.load_by_name(name)

        # if neither a canonical or alias exists, create the records. Remember to check for
        # canonical / alias status.
        else:
            lp = LigandProperities()
            lp.save()
            l = Ligand()
            l.properities = lp
            l.name = str(name)
            l.canonical = True
            l.ambigious_alias = False
            try:
                l.save()
                l.load_by_name(str(name))
            except IntegrityError:
                l = Ligand.objects.get(name=str(name), canonical=True)
    else:
        l = None

    return l

#def fetch_chembl_refs(lig_chembl_id, target_accesion):

#    target_id = new_client.target.filter(accession=target_accesion)

#    assay = new_client.assay.filter(target_id=target_id, compound=lig_chembl_id)

#    refs = [x['document_chembl_id'] for x in assay]
#    #https://www.ebi.ac.uk/chembl/doc/inspect/CHEMBL2766014

#    return refs
