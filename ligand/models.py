from django.db import models

from common.models import WebResource
from common.models import WebLink

from urllib.request import urlopen
import json


class Ligand(models.Model):
    ligand_type = models.ForeignKey('LigandType', null=True)
    web_links = models.ManyToManyField('common.WebLink')
    name = models.TextField()
    smiles = models.TextField(null=True)
    inchikey = models.CharField(max_length=50, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'ligand'


    def load_by_name(self, name):
        # fetch ligand info from pubchem
        pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + name + '/json'
        
        try:
            req = urlopen(pubchem_url)
            pubchem = json.loads(req.read().decode('UTF-8'))
        except:
            return

        # weblink
        pubchem_id = pubchem['PC_Compounds'][0]['id']['id']['cid']
        try:
            web_resource = WebResource.objects.get(slug='pubchem')
        except:
            # abort if pdb resource is not found
            raise Exception('PubChem resource not found, aborting!')
        wl, created = WebLink.objects.get_or_create(index=pubchem_id, web_resource=web_resource)
        self.web_links.add(wl)

        # ligand type
        self.ligand_type, created = LigandType.objects.get_or_create(slug='sm', defaults={'name':'Small molecule'})

        # SMILES
        for prop in pubchem['PC_Compounds'][0]['props']:
            if prop['urn']['label'] == 'SMILES' and prop['urn']['name'] == 'Canonical':
                self.smiles = prop['value']['sval']
                break

        # InChIKey
        for prop in pubchem['PC_Compounds'][0]['props']:
            if prop['urn']['label'] == 'InChIKey':
                self.inchikey = prop['value']['sval']
                break


class LigandType(models.Model):
    slug = models.SlugField(max_length=50)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_type'


class LigandAlias(models.Model):
    ligand = models.ForeignKey('Ligand')
    name = models.TextField()

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_alias'


class LigandRole(models.Model):
    slug = models.SlugField(max_length=50)
    name = models.CharField(max_length=100)
    
    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_role'