from django.db import models
from django.utils.text import slugify
from django.db import IntegrityError

from common.models import WebResource
from common.models import WebLink
from common.tools import fetch_from_web_api

from urllib.request import urlopen, quote
import json
import yaml
import logging

class Ligand(models.Model):
    properities = models.ForeignKey('LigandProperities', on_delete=models.CASCADE)
    name = models.TextField()
    canonical = models.NullBooleanField()
    ambigious_alias = models.NullBooleanField() #required to flag 'safe' alias, eg one parent

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand'
        unique_together = ('name', 'canonical')

    def load_by_gtop_id(self, ligand_name, gtop_id, ligand_type):
        logger = logging.getLogger('build')

        # get the data from cache or web services
        cache_dir = ['guidetopharmacology', 'ligands']
        url = 'http://www.guidetopharmacology.org/services/ligands/$index'
        gtop = fetch_from_web_api(url, gtop_id, cache_dir)
        
        if gtop:
            # get name from response
            ligand_name = gtop['name']
            if ligand_name=='11-<i>cis</i>-retinal':
                ligand_name = 'retinal'

        # does a ligand by this name already exists?
        try:
            existing_ligand = Ligand.objects.get(name=ligand_name, canonical=True)
            return existing_ligand
        except Ligand.DoesNotExist:
            web_resource = False
            
            if gtop_id:
                # gtoplig webresource
                web_resource = WebResource.objects.get(slug='gtoplig')
            
            return self.update_ligand(ligand_name, {}, ligand_type, web_resource, gtop_id)

    def load_from_pubchem(self, lookup_type, pubchem_id, ligand_type, ligand_title=False):
        logger = logging.getLogger('build')

        # if ligand title is specified, use that as the name
        if ligand_title:
            ligand_name = ligand_title

        # otherwise, fetch ligand name from pubchem
        else:
            # check cache
            cache_dir = ['pubchem', 'cid', 'synonyms']
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/$index/synonyms/json'.format(lookup_type)
            pubchem = fetch_from_web_api(url, pubchem_id, cache_dir)
            ##print (pubchem)
            
            # get name from response
            try:
                ligand_name = pubchem['InformationList']['Information'][0]['Synonym'][0]
            except:
                ## Some compounds do not have a name but are still a valid pubchem entry. (Peptides)
                logger.warning('Ligand {} does not have a name in PubChem'.format(pubchem_id))
                ligand_name = lookup_type + ' ' + pubchem_id
                # return None

        # fetch ligand properties from pubchem
        properties = {}
        
        # check cache
        cache_dir = ['pubchem', 'cid', 'property']
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/$index/property/CanonicalSMILES,InChIKey,MolecularWeight,HBondDonorCount,HBondAcceptorCount,XLogP,RotatableBondCount/json'.format(lookup_type)
        pubchem = fetch_from_web_api(url, pubchem_id, cache_dir)
        # get properties from response
        if pubchem==False:
            logger.warning('Ligand {} not found in PubChem'.format(pubchem_id))
            return None

        if pubchem['PropertyTable']['Properties'][0]:   
            if 'HBondAcceptorCount' in pubchem['PropertyTable']['Properties'][0] :
                properties['hacc'] =  pubchem['PropertyTable']['Properties'][0]['HBondAcceptorCount']
            if 'HBondDonorCount' in pubchem['PropertyTable']['Properties'][0] :
                properties['hdon'] =  pubchem['PropertyTable']['Properties'][0]['HBondDonorCount']
            if 'XLogP' in pubchem['PropertyTable']['Properties'][0] :
                properties['logp'] =  pubchem['PropertyTable']['Properties'][0]['XLogP']
            if 'RotatableBondCount' in pubchem['PropertyTable']['Properties'][0] :
                properties['rotatable_bonds'] =  pubchem['PropertyTable']['Properties'][0]['RotatableBondCount']
            if 'MolecularWeight' in pubchem['PropertyTable']['Properties'][0] :
                properties['mw'] = pubchem['PropertyTable']['Properties'][0]['MolecularWeight']
        try:
            
            properties['smiles'] =  pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
            properties['inchikey'] =  pubchem['PropertyTable']['Properties'][0]['InChIKey']

        except:
            logger.warning('Ligand {} not found in PubChem'.format(pubchem_id))
            return None

        # pubchem webresource
        web_resource = WebResource.objects.get(slug='pubchem')
        #print (web_resource)
        
        # does a ligand with this canonical name already exist
        try:
            return Ligand.objects.get(name=ligand_name, canonical=True)
            # FIXME check inchikey
        except Ligand.DoesNotExist:
            pass # continue

        # does a (canonical) ligand with this inchikey already exist?
        try:
            existing_lp = LigandProperities.objects.get(inchikey=properties['inchikey'])
            self.properities = existing_lp
            self.name = ligand_name
            self.canonical = False
            self.ambigious_alias = False
            
            try:
                self.save()
                return self
            except IntegrityError:
                return Ligand.objects.get(name=ligand_name, canonical=False)
        except LigandProperities.DoesNotExist:
            return self.update_ligand(ligand_name, properties, ligand_type, web_resource, pubchem_id)

    def update_ligand(self, ligand_name, properties, ligand_type, web_resource=False, web_resource_index=False):
        lp = LigandProperities.objects.create(ligand_type=ligand_type)

        # assign properties
        for prop in properties:
            setattr(lp, prop, properties[prop])

        # assign web link
        if web_resource and web_resource_index:
            try:
                wl, created = WebLink.objects.get_or_create(index=web_resource_index, web_resource=web_resource)
            except IntegrityError:
                wl = Weblink.objects.get(index=web_resource_index, web_resource=web_resource)
            lp.web_links.add(wl)

        # try saving the properties, catch IntegrityErrors due to concurrent processing
        try:
            lp.save()
        except IntegrityError:
            lp = LigandProperities.objects.get(inchikey=properties['inchikey'])

        self.name = ligand_name
        self.canonical = True
        self.ambigious_alias = False
        self.properities = lp
        
        try:
            self.save()
            return self
        except IntegrityError:
            return Ligand.objects.get(name=ligand_name, canonical=True)

    def load_by_name(self, name):
        logger = logging.getLogger('build')
        # fetch ligand info from pubchem - start by getting name and 'canonical' name
        pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + name + '/synonyms/TXT'
        if self.properities.inchikey: #if inchikey has been added use this -- especially to avoid updating a wrong inchikey to a synonym. 
            pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/InchiKey/' + self.properities.inchikey + '/synonyms/TXT'
        try:
            req = urlopen(pubchem_url)
            pubchem = req.read().decode('UTF-8').splitlines()
            pubchem_name = pubchem[0]
        except: #name not matched in pubchem 
            if self.properities.inchikey: #if inchikey has been added for check this
                pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/InchiKey/' + self.properities.inchikey + '/synonyms/TXT'
                try:
                    req = urlopen(pubchem_url)
                    pubchem = req.read().decode('UTF-8').splitlines()
                    pubchem_name = pubchem[0]
                except: #name not matched in pubchem - exit cos something is wrong
                    logger.info('Ligand not found by InchiKey in pubchem: ' + str(self.properities.inchikey))
                    return
            else: #if name not found and no inchikey, then no point in looking further
                logger.info('Ligand not found in pubchem by name (Consider renaming): ' + str(name))
                return

        pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + quote(pubchem_name) + '/property/CanonicalSMILES,InChIKey/json'

        if self.properities.inchikey:
            pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inChiKey/' + self.properities.inchikey + '/property/CanonicalSMILES,InChIKey/json'

        try:
            req = urlopen(pubchem_url)
            pubchem = json.loads(req.read().decode('UTF-8'))
        except: #JSON failed
            return

        # weblink
        pubchem_id = pubchem['PropertyTable']['Properties'][0]['CID']
        try:
            web_resource = WebResource.objects.get(slug='pubchem')
        except:
            # abort if pdb resource is not found
            raise Exception('PubChem resource not found, aborting!')
        
        pubchem_inchikey = ''
        pubchem_smiles = ''

        # SMILES
        pubchem_smiles = pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']

        # InChIKey
        pubchem_inchikey = pubchem['PropertyTable']['Properties'][0]['InChIKey']

        try: #now that we have inchikey, try and see if it exists in DB
            existing_lp = LigandProperities.objects.get(inchikey=pubchem_inchikey)
            self.properities = existing_lp
        except:
            wl, created = WebLink.objects.get_or_create(index=pubchem_id, web_resource=web_resource)
            self.properities.web_links.add(wl)
            # ligand type
            self.properities.ligand_type, created = LigandType.objects.get_or_create(slug='sm', defaults={'name':'Small molecule'})
            self.properities.inchikey = pubchem_inchikey
            self.properities.smiles = pubchem_smiles
            self.properities.save()

        if pubchem_name.lower()!=name.lower(): #if not canonical name
            logger.info("Updating canonical flag to Pubchem. PubChem canonical: "+pubchem_name +". DB canonical: "+ name)
            self.canonical = False
            try:
                self.save()
            except IntegrityError:
                logger.error("FAILED SAVING LIGAND, duplicate?")
            canonical_entry = Ligand.objects.filter(name=pubchem_name, properities__inchikey=pubchem_inchikey) #NEED TO CHECK BY INCHI KEY - SOME CANONICAL NAMES HAVE MANY ICHIKEYS (DOXEPIN)
            if canonical_entry.exists():
                return
            else: #insert the 'canonical' entry
                try:
                    canonical_entry = Ligand()
                    canonical_entry.name = pubchem_name
                    canonical_entry.canonical = True
                    canonical_entry.properities = self.properities
                    canonical_entry.save()
                except IntegrityError:
                    logger.error("FAILED SAVING CANONICAL LIGAND, duplicate? "+pubchem_name+" "+name)
                    print("FAILED SAVING CANONICAL LIGAND, duplicate? "+pubchem_name+" "+name)


class LigandProperities(models.Model):
    ligand_type = models.ForeignKey('LigandType', null=True, on_delete=models.CASCADE)
    web_links = models.ManyToManyField('common.WebLink')
    #vendor_links = models.ManyToManyField('common.WebLink', related_name='vendors')
    smiles = models.TextField(null=True)
    inchikey = models.CharField(max_length=50, null=True, unique=True)
    #vendors = models.ManyToManyField('LigandVenderLink')

    mw = models.DecimalField(max_digits=15, decimal_places=3, null=True)
    rotatable_bonds =  models.SmallIntegerField(null=True)
    hacc =  models.SmallIntegerField( null=True)
    hdon =  models.SmallIntegerField( null=True)
    logp = models.DecimalField(max_digits=10, decimal_places=3, null=True)
    

    def __str__(self):
        return self.inchikey

    class Meta():
        db_table = 'ligand_properities'


    
class LigandType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_type'


class LigandRole(models.Model):
    slug = models.SlugField(max_length=50, unique=True)
    name = models.CharField(max_length=100)
    
    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_role'
        
        
class ChemblAssay(models.Model):
     #slug = models.SlugField(max_length=50, unique=True)
    web_links = models.ManyToManyField('common.WebLink')
    assay_id = models.CharField(max_length=50, unique = True)
    
    
        
    def __str__(self):
        return self.assay_id
    
    class Meta():
        db_table = 'chembl_assays'
        
        
class AssayExperiment(models.Model):
    
    ligand = models.ForeignKey('Ligand', on_delete=models.CASCADE)
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    assay = models.ForeignKey('ChemblAssay', on_delete=models.CASCADE)
    assay_type = models.CharField(max_length=10)
    assay_description = models.TextField(max_length=1000)
    pchembl_value = models.DecimalField(max_digits=9, decimal_places=3)

    published_value = models.DecimalField(max_digits=9, decimal_places=3)
    published_relation = models.CharField(max_length=10)
    published_type = models.CharField(max_length=20)
    published_units = models.CharField(max_length=20)

    standard_value = models.DecimalField(max_digits=9, decimal_places=3)
    standard_relation = models.CharField(max_length=10)
    standard_type = models.CharField(max_length=20)
    standard_units = models.CharField(max_length=20)
    

    
    class Meta():
        unique_together = ('ligand', 'protein', 'assay')
    
    

class LigandVendors(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    name = models.CharField(max_length=200, default='')
    url = models.TextField(null=True)

class LigandVendorLink(models.Model):
    vendor = models.ForeignKey('LigandVendors', on_delete=models.CASCADE)
    lp = models.ForeignKey('LigandProperities', related_name='vendors', on_delete=models.CASCADE)
    url = models.CharField(max_length=300)  #SourceRecordURL
    vendor_external_id = models.CharField(max_length=300) #RegistryID
    sid = models.CharField(max_length=200, unique=True) #SID
    
