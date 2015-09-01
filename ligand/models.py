from django.db import models

from common.models import WebResource
from common.models import WebLink

from urllib.request import urlopen, quote
import json
import logging

class Ligand(models.Model):
    properities = models.ForeignKey('LigandProperities')
    name = models.TextField()
    canonical = models.NullBooleanField()
    ambigious_alias = models.NullBooleanField() #required to flag 'safe' alias, eg one parent 

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand'

    def update_by_PubChemId(self, pubchem_id):
        #IUPACName,
        pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' + pubchem_id + '/property/CanonicalSMILES,InChIKey/json'

        #pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + quote(pubchem_name) + '/json'
        try:
            req = urlopen(pubchem_url)
            pubchem = json.loads(req.read().decode('UTF-8'))
        except: #JSON failed
            return

        pubchem_smiles = pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
        pubchem_inchikey = pubchem['PropertyTable']['Properties'][0]['InChIKey']



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

        wl, created = WebLink.objects.get_or_create(index=pubchem_id, web_resource=web_resource)
        self.properities.web_links.add(wl)
        # ligand type
        self.properities.ligand_type, created = LigandType.objects.get_or_create(slug='sm', defaults={'name':'Small molecule'})
        self.properities.inchikey = pubchem_inchikey
        self.properities.smiles = pubchem_smiles
        self.properities.save()

        if pubchem_name!=name: #if not canonical name
            logger.info("Updating canonical flag to Pubchem. PubChem canonical: "+pubchem_name +". DB canonical: "+ name)
            self.canonical = False
            self.save()
            canonical_entry = Ligand.objects.filter(name=pubchem_name, properities__inchikey=pubchem_inchikey) #NEED TO CHECK BY INCHI KEY - SOME CANONICAL NAMES HAVE MANY ICHIKEYS (DOXEPIN)
            if canonical_entry.exists():
                return
            else: #insert the 'canonical' entry
                canonical_entry = Ligand()
                canonical_entry.name = pubchem_name
                canonical_entry.canonical = True
                canonical_entry.properities = self.properities
                canonical_entry.save()


class LigandProperities(models.Model):
    ligand_type = models.ForeignKey('LigandType', null=True)
    web_links = models.ManyToManyField('common.WebLink')
    smiles = models.TextField(null=True)
    inchikey = models.CharField(max_length=50, null=True, unique=True)


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
    slug = models.SlugField(max_length=50)
    name = models.CharField(max_length=100)
    
    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_role'