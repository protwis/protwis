from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from common.tools import test_model_updates
from build.management.commands.parse_excel_annotations import Command as ParseExcel
from common.models import Publication, PublicationJournal, WebLink, WebResource, Citation, DefaultCitation

import django.apps
import datetime
import logging
import os
import yaml
import pprint
import re
from collections import OrderedDict



class Command(ParseExcel):
    help = 'Basic functions for build scripts'

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            help='Purge citation entries')
    REF_HEADER = 'Reference DOI'
    re_ref_header = re.compile(r'^'+re.escape(REF_HEADER)+r'(\s+(\d+))?$',flags=re.IGNORECASE)
    source_file = os.sep.join([settings.DATA_DIR, 'common_data', 'Site_References.xlsx'])
    references_yaml = os.sep.join([settings.DATA_DIR, 'common_data', 'references.yaml'])
    default_references_yaml = os.sep.join([settings.DATA_DIR, 'common_data', 'default_references.yaml'])
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
        with open(self.default_references_yaml, 'r') as refs_yaml:
            default_refs = yaml.load(refs_yaml, Loader=yaml.FullLoader)
        wr = WebResource.objects.get(slug='doi')
        pubjournal = None

        default_refs_pub_dict = {}
        for default_ref_name, vals in default_refs.items():
            pubs = []
            for ypub in sorted(vals, key = lambda x:x['Order']):
                doi = ypub['DOI']
                pub = False
                print('doi:',ypub['DOI'])

                title = ypub['Title']
                there_is_title = len(title) > 0
                year = ypub['Year']
                created_pj = False
                pubjournal = None
                if there_is_title:
                    pubjournal, created_pj = PublicationJournal.objects.get_or_create(defaults={"name": ypub["Journal"], 'slug': slugify(ypub['Journal'])}, name__iexact=ypub["Journal"])
                pub = self.create_publication(doi, wr, pubjournal,force=there_is_title,delete_pj=created_pj)
                edit = False
                if there_is_title and not pub.title:
                    pub.title = title
                    edit = True
                if year and not pub.year:
                    pub.year = year
                    edit = True
                if edit:
                    pub.save()

                print(pub)
                if pub:
                    pubs.append(pub)
            default_refs_pub_dict[default_ref_name] = pubs
            dcit, created = DefaultCitation.objects.get_or_create(main=default_ref_name)
            for pub in pubs:
                dcit.publication.add(pub)
            dcit.save()
            self.logger.info('Created default citation: {}'.format(dcit))

        for url, vals in refs.items():
            if vals['Default']!=0:
                main = vals['Default']
            else:
                main = None
            page = vals['Page']

            if vals['Video']!='':
                vid = vals['Video']
            else:
                vid = None
            if main is not None:
                default = main
            else:
                default = vals['Menu']
            pubs = []
            for ypub in vals['Publication']:
                doi = ypub['DOI']
                title = ypub['Title']
                year = ypub['Year']
                there_is_title = len(title) > 0
                pub = False
                print('juas juas1')
                
                print('juas juas2')

                created_pj = False
                pubjournal = None
                if there_is_title:
                    pubjournal, created_pj = PublicationJournal.objects.get_or_create(defaults={"name": ypub["Journal"], 'slug': slugify(ypub['Journal'])}, name__iexact=ypub["Journal"])

                pub = self.create_publication(doi, wr, pubjournal,force=there_is_title,delete_pj=created_pj)

                edit = False
                if there_is_title and not pub.title:
                    pub.title = title
                    edit = True
                if year and not pub.year:
                    pub.year = year
                    edit = True
                if edit:
                    pub.save()
                if pub:
                    pubs.append(pub)
            print('juas juas4')
            cit, created = Citation.objects.get_or_create(defaults={'video':vid,'main':main,'docs':None,'page_name':page},url=url)

            for pub in pubs:
                print('juas pub', pub)
                cit.publication.add(pub)
            print('juas juas5')
            cit.save()
            self.logger.info('Created citation: {}'.format(cit))

    def purge_citations(self):
        Citation.objects.all().delete()
        DefaultCitation.objects.all().delete()

    def write_yaml_refs(self, data):
        '''Write parsed site references data to yaml file'''
        
        re_ref_header = self.re_ref_header
        cit_dict = OrderedDict()
        error_ref_header = None
        for i, vals in data['References'].items():
            url = vals['URL']
            if url == 'broken':
                continue


            nref_set = set()
            ref_header_dict = {}
            ref_header = len(self.REF_HEADER)
            for k in vals.keys():
                m = re_ref_header.match(k.strip())
                if m:
                    num = m.group(2)
                    if num is None:
                        num = 1                  
                    num = int(num)
                    if num in nref_set:
                        error_num = num
                        error_ref_header = k
                    else:
                        nref_set.add(num)
                        ref_header_dict[num] = k

            if error_ref_header is not None:
                raise ValueError('Invalid DOI reference header: '+k)
            
            publication_list = []
            for num,h in sorted(ref_header_dict.items(),key = lambda x:x[0]):
                if num == 1:
                    publication_list.append({'DOI':self.parse_DOI(vals[h]),
                                            'Journal':vals['Journal'],
                                            'Title':vals['Title'],
                                            'Year':vals['Year']})
                else:
                    publication_list.append({'DOI':self.parse_DOI(vals[h]),
                        'Journal':vals['Journal'+' '+str(num)],
                        'Title':vals['Title'+' '+str(num)],
                        'Year':vals['Year'+' '+str(num)]})

            publication_list = [p for p in publication_list 
                                if not ((p['DOI'] == '' or p['DOI'] is None) and
                                       (p['Journal'] == '' or p['Journal'] is None) == '' and
                                       (p['Title'] == '' or p['Title'] is None) and
                                        p['Year'] == '' or p['Year'] is None)]
            

            cit_dict[url] = {'Menu':vals['Menu'], 'Section':vals['Section'], 'Menu Section':vals['Menu Section'],
                             'Menu subsection':vals['Menu subsection'], 'Page':vals['Page'], 'Video':vals['Video URL'],
                             'Default':vals['Default reference'],
                             'Publication':publication_list}                      

        with open(self.references_yaml, 'w') as f:
            yaml.dump(cit_dict, f, indent=4)

        cit_dict = OrderedDict()
        for i, vals in data['DefaultReferences'].items():
            default_reference = vals['Default reference']
            if default_reference not in cit_dict:
                cit_dict[default_reference] = []
            cit_dict[default_reference].append({'DOI':self.parse_DOI(vals[self.REF_HEADER]),
                                            'Journal':vals['Journal'],
                                            'Title':vals['Title'],
                                            'Year':vals['Year'],
                                            'Order':int(vals['Order'])})
        
        with open(self.default_references_yaml, 'w') as f:
            yaml.dump(cit_dict, f, indent=4)

        

    def parse_DOI(self, url):
        if 'https://doi.org/' in url:
            return url.split('https://doi.org/')[1].upper()
        elif 'http://doi.org/' in url:
            return url.split('http://doi.org/')[1].upper()
        elif 'https://dx.doi.org/' in url:
            return url.split('https://dx.doi.org/')[1].upper()
        elif 'http://dx.doi.org/' in url:
            return url.split('http://dx.doi.org/')[1].upper()
        else:
            return url.upper()

    def create_publication(self, doi, wr, pubjournal,force=False,delete_pj=False):
        '''Create WebLink and Publication objects'''
        if doi!='':
            try:
                pub = Publication.objects.get(web_link__index=doi, web_link__web_resource=wr)
            except Publication.DoesNotExist as e:
                pub = Publication.get_or_create_from_doi(doi)

                if not pub.journal and pubjournal:
                    pub.journal = pubjournal
                    pub.save()
                    delete_pj = False

                self.logger.info('Created Publication:'+str(pub))
            if not pub.journal and pubjournal:
                pub.journal = pubjournal
                pub.save()
                delete_pj = False
            if delete_pj:
                pubjournal.delete()
            return pub
        elif force:
            pub, create = Publication.objects.get_or_create(journal=pubjournal)
            return pub
        else:
            return False
