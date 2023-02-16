from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from news.models import News
from common.tools import test_model_updates
import logging
import os
import yaml
import django.apps

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)

    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    news_data_dir = os.sep.join([settings.DATA_DIR, 'news'])

    def handle(self, *args, **options):
        functions = [
            'create_news',
        ]

        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        test_model_updates(self.all_models, self.tracker, check=True)

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
                    ds = yaml.load(f, Loader=yaml.FullLoader)
                    doc, created = News.objects.get_or_create(image=ds['image'], date=ds['date'])
                    if created:
                        self.logger.info('Created news for ' + ds['date'])

                with open(source_file_path[:-4]+'html', 'r') as h:
                    doc.html = h.read()
                    doc.save()
                    if created:
                        self.logger.info('Created html news for ' + ds['date'])

        self.logger.info('COMPLETED CREATING NEWS')
