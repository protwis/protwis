from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from common.tools import test_model_updates
from build.management.commands.parse_excel_annotations import Command as ParseExcel
from common.models import Publication, PublicationJournal, WebLink, WebResource, Citation

import django.apps
import datetime
import logging
import os
import yaml
import pprint
from collections import OrderedDict

class Command(ParseExcel):
    help = 'Basic functions for build scripts'

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            help='Purge citation entries')

    source_file = os.sep.join([settings.DATA_DIR, 'common_data', 'Site_References.xlsx'])
    references_yaml = os.sep.join([settings.DATA_DIR, 'common_data', 'references.yaml'])
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def handle(self, *args, **options):
        if options['purge']:
            self.purge_citations()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)
        self.logger.info('CREATING CITATIONS')
        self.logger.info('Parsing file ' + self.source_file)
        data = self.parse_excel(self.source_file)
        self.logger.info('Writing data to ' + self.references_yaml)
        self.write_yaml_refs(data)
        self.logger.info('Parsing {} and creating citation entries'.format(self.references_yaml))
        self.create_citations()
        test_model_updates(self.all_models, self.tracker, check=True)

    def create_citations(self):
        '''Create citation objects'''
        with open(self.references_yaml, 'r') as refs_yaml:
            refs = yaml.load(refs_yaml, Loader=yaml.FullLoader)
        wr = WebResource.objects.get(slug='doi')
        pubjournal = None

        # Create main publication first for empty publication cells for non-published tools
        for url, vals in refs.items():
            if vals['Default']=='GPCRdb':
                main_gpcrdb_pub = Publication.get_or_create_from_doi(vals['DOI'])
            if vals['Default']=='GproteinDb':
                main_gproteindb_pub = Publication.get_or_create_from_doi(vals['DOI'])

        for url, vals in refs.items():
            doi = vals['DOI']
            pub = False
            if vals['Journal'] in ['Preprint at Research Square', 'Submitted']:
                pubjournal, created = PublicationJournal.objects.get_or_create(defaults={"name": vals["Journal"], 'slug': slugify(vals['Journal'])}, name__iexact=vals["Journal"])
                pub = self.create_publication(doi, wr, pubjournal)
            elif len(doi) > 0:
                pub = Publication.get_or_create_from_doi(vals['DOI'])

            page = vals['Page']
            if not pub:
                if vals['Menu']=='Gproteindb':
                    pub = main_gproteindb_pub
                else:
                    pub = main_gpcrdb_pub

            if vals['Video']!='':
                vid = vals['Video']
            else:
                vid = None

            if vals['Default']!=0:
                main = vals['Default']
            else:
                main = None

            cit, created = Citation.objects.get_or_create(publication=pub, url=url, video=vid, main=main, docs=None, page_name=page)
            cit.save()
            self.logger.info('Created citation: {}'.format(cit))

    def purge_citations(self):
        Citation.objects.all().delete()

    def write_yaml_refs(self, data):
        '''Write parsed site references data to yaml file'''
        cit_dict = OrderedDict()
        for i, vals in data['References'].items():
            cit_dict[vals['URL']] = {'Menu':vals['Menu'], 'Section':vals['Section'], 'Page':vals['Page'], 'Video':vals['Video URL'],
                                     'DOI':self.parse_DOI(vals['Reference DOI']), 'Default':vals['Default reference'], 'Journal':vals['Journal'],
                                     'Title':vals['Title'], 'Year':vals['Year']}
        with open(self.references_yaml, 'w') as f:
            yaml.dump(cit_dict, f, indent=4)

    def parse_DOI(self, url):
        if 'https://doi.org/' in url:
            return url.split('https://doi.org/')[1].upper()
        else:
            return url.upper()

    def create_publication(self, doi, wr, pubjournal):
        '''Create WebLink and Publication objects'''
        if doi!='':
            try:
                pub = Publication.objects.get(web_link__index=doi, web_link__web_resource=wr)
            except Publication.DoesNotExist as e:
                pub = Publication.get_or_create_from_doi(doi)

                if not pub.journal and pubjournal:
                    pub.journal = pubjournal
                    pub.save()
                self.logger.info('Created Publication:'+str(pub))
            return pub
        else:
            return None
