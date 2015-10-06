from django.db import models
from django.utils.text import slugify
from django.db import IntegrityError

from common.models import WebResource
from common.models import WebLink
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api

from urllib.request import urlopen, quote
import json
import yaml
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
        unique_together = ('name', 'canonical')

    def load_by_gtop_id(self, ligand_name, gtop_id, ligand_type):
        logger = logging.getLogger('build')

        # check whether this data is cached
        cache_dir = ['guidetopharmacology', 'ligands']
        gtop = fetch_from_cache(cache_dir, str(gtop_id))

        if gtop:
            logger.info('Fetched {}/{} from cache'.format('/'.join(cache_dir), gtop_id))
        else:
            logger.info('No cached entry for {}/{}'.format('/'.join(cache_dir), gtop_id))
            
            # fetch synomyms
            gtop_url = 'http://www.guidetopharmacology.org/services/ligands/' + str(gtop_id)

            try:
                req = fetch_from_web_api(gtop_url)
                if req:
                    gtop = json.loads(req.read().decode('UTF-8'))

                    # save to cache
                    save_to_cache(cache_dir, str(gtop_id), gtop)
                    logger.info('Saved entry for {}/{} in cache'.format('/'.join(cache_dir), gtop_id))
            except:
                logger.error('Failed fetching properties of ligand with GuideToPharmacology ID {}'.format(gtop_id))
        
        if gtop:
            # get name from response
            ligand_name = gtop['name']

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

    def load_by_pubchem_id(self, pubchem_id, ligand_type, ligand_title):
        logger = logging.getLogger('build')

        # if ligand title is specified, use that as the name
        if ligand_title:
            ligand_name = ligand_title

        # otherwise, fetch ligand name from pubchem
        else:
            # check cache
            cache_dir = ['pubchem', 'cid', 'synonyms']
            pubchem = fetch_from_cache(cache_dir, str(pubchem_id))

            if pubchem:
                logger.info('Fetched {}/{} from cache'.format('/'.join(cache_dir), pubchem_id))
            else:
                logger.info('No cached entry for {}/{}'.format('/'.join(cache_dir), pubchem_id))

                pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' + str(pubchem_id) \
                    + '/synonyms/json'
                try:
                    req = fetch_from_web_api(pubchem_url)
                    if req:
                        pubchem = json.loads(req.read().decode('UTF-8'))

                        # save to cache
                        save_to_cache(cache_dir, str(pubchem_id), pubchem)
                        logger.info('Saved entry for {}/{} in cache'.format('/'.join(cache_dir), pubchem_id))
                except:
                    logger.error('Error fetching ligand {} from PubChem'.format(pubchem_id))
                    return
            
            # get name from response
            ligand_name = pubchem['InformationList']['Information'][0]['Synonym'][0]

        # fetch ligand properties from pubchem
        properties = {}
        
        # check cache
        cache_dir = ['pubchem', 'cid', 'property']
        pubchem = fetch_from_cache(cache_dir, str(pubchem_id))

        if pubchem:
            logger.info('Fetched {}/{} from cache'.format('/'.join(cache_dir), pubchem_id))
        else:
            logger.info('No cached entry for {}/{}'.format('/'.join(cache_dir), pubchem_id))
            
            pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' + str(pubchem_id) \
                + '/property/CanonicalSMILES,InChIKey/json'
            try:
                req = fetch_from_web_api(pubchem_url)
                if req:
                    pubchem = json.loads(req.read().decode('UTF-8'))

                    # save to cache
                    save_to_cache(cache_dir, str(pubchem_id), pubchem)
                    logger.info('Saved entry for {}/{} in cache'.format('/'.join(cache_dir), pubchem_id))
            except:
                logger.error('Error fetching ligand {} from PubChem'.format(pubchem_id))
                return
        
        # get properties from reponse
        properties['smiles'] =  pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
        properties['inchikey'] =  pubchem['PropertyTable']['Properties'][0]['InChIKey']

        # pubchem webresource
        web_resource = WebResource.objects.get(slug='pubchem')

        # does a ligand with this inchikey already exists?
        try:
            existing_ligand = Ligand.objects.get(properities__inchikey=properties['inchikey'], canonical=True)
            self.properities = existing_ligand.properities
            self.name = ligand_name
            self.canonical = False
            self.ambigious_alias = False
            
            try:
                self.save()
                return self
            except IntegrityError:
                return Ligand.objects.get(name=ligand_name, canonical=False)
        except Ligand.DoesNotExist:
            return self.update_ligand(ligand_name, properties, ligand_type, web_resource, pubchem_id)

    def update_ligand(self, ligand_name, properties, ligand_type, web_resource=False, web_resource_index=False):
        lp = LigandProperities()
        lp.ligand_type = ligand_type
        lp.save()

        # assign properties
        for prop in properties:
            setattr(lp, prop, properties[prop])

        # assign web link
        if web_resource and web_resource_index:
            wl, created = WebLink.objects.get_or_create(index=web_resource_index, web_resource=web_resource)
            lp.web_links.add(wl)

        lp.save()
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