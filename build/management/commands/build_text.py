from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from documentation.models import Documentation
from news.models import News
from pages.models import Pages

import logging
import os
import yaml

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)

    documentation_data_dir = os.sep.join([settings.DATA_DIR, 'documentation'])
    news_data_dir = os.sep.join([settings.DATA_DIR, 'news'])
    pages_data_dir = os.sep.join([settings.DATA_DIR, 'pages'])

    def handle(self, *args, **options):
        functions = [
            'create_documentation',
            'create_news',
            'create_pages',
        ]

        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def create_documentation(self):
        self.logger.info('CREATING DOCUMENTATION')

        # clear existing documentation
        Documentation.objects.all().delete()
        
        # what files should be parsed?
        filenames = os.listdir(self.documentation_data_dir)

        for source_file in filenames:
            self.logger.info('Parsing file ' + source_file)
            source_file_path = os.sep.join([self.documentation_data_dir, source_file])

            if source_file.endswith('.yaml'):
                with open(source_file_path, 'r') as f:
                    ds = yaml.load(f)
                    doc, created = Documentation.objects.get_or_create(title=ds['title'],
                        description=ds['description'], image=ds['image'] ,)
                    if created:
                        self.logger.info('Created documentation for ' + ds['title'])
                
                with open(source_file_path[:-4]+'html', 'r') as h:
                    doc.html = h.read()
                    doc.save()
                    if created:
                        self.logger.info('Created html documentation for ' + ds['title'])

        self.logger.info('COMPLETED CREATING DOCUMENTATION')

    def create_pages(self):
        self.logger.info('CREATING PAGES')

        # clear existing pages
        Pages.objects.all().delete()
        
        # what files should be parsed?
        filenames = sorted(os.listdir(self.pages_data_dir),key=str.lower)

        for source_file in filenames:
            self.logger.info('Parsing file ' + source_file)
            source_file_path = os.sep.join([self.pages_data_dir, source_file])

            if source_file.endswith('.yaml'):
                with open(source_file_path, 'r') as f:
                    ds = yaml.load(f)
                    doc, created = Pages.objects.get_or_create(title=ds['title'])
                    if created:
                        self.logger.info('Created pages for ' + ds['title'])

                with open(source_file_path[:-4]+'html', 'r') as h:
                    doc.html = h.read()
                    doc.save()
                    if created:
                        self.logger.info('Created html pages for ' + ds['title'])

        self.logger.info('COMPLETED CREATING PAGES')

    def create_news(self):
        self.logger.info('CREATING NEWS')

        # clear existing news
        News.objects.all().delete()
        
        # what files should be parsed?
        filenames = os.listdir(self.news_data_dir)

        for source_file in filenames:
            self.logger.info('Parsing file ' + source_file)
            source_file_path = os.sep.join([self.news_data_dir, source_file])

            if source_file.endswith('.yaml'):
                with open(source_file_path, 'r') as f:
                    ds = yaml.load(f)
                    doc, created = News.objects.get_or_create(image=ds['image'], date=ds['date'])
                    if created:
                        self.logger.info('Created news for ' + ds['date'])
                
                with open(source_file_path[:-4]+'html', 'r') as h:
                    doc.html = h.read()
                    doc.save()
                    if created:
                        self.logger.info('Created html news for ' + ds['date'])

        self.logger.info('COMPLETED CREATING NEWS')
