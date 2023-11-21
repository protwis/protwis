from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from common.models import WebResource, WebLink
from protein.models import Protein
from optparse import make_option
from common.tools import test_model_updates
import django.apps
import logging
import shlex
import os

class Command(BaseCommand):
    help = 'Reads source data and creates links to other databases'

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

    # source file directory
    links_data_dir = os.sep.join([settings.DATA_DIR, 'protein_data', 'links'])
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        try:
            self.create_links(filenames)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def create_links(self, filenames=False):
        self.logger.info('CREATING LINKS')

        # read source files
        if not filenames:
            filenames = os.listdir(self.links_data_dir)

        for filename in filenames:
            # get the resource name from the filename (files should be named "resource_slug.txt")
            resource_slug = filename.split('.')[0]

            # find resource
            try:
                resource = WebResource.objects.get(slug=resource_slug)
            except WebResource.DoesNotExist:
                self.logger.error('Web resource {} does not exits'.format(resource_slug))
                continue

            # counter for created links
            num_created_links = 0

            filepath = os.sep.join([self.links_data_dir, filename])
            with open(filepath, "r", encoding='UTF-8') as f:
                self.logger.info('Parsing file {}'.format(filepath))

                for row in f:
                    split_row = shlex.split(row)
                    accession = split_row[0]
                    link_index = split_row[1]

                    # fetch protein
                    try:
                        p = Protein.objects.get(accession=accession)
                    except Protein.DoesNotExist:
                        self.logger.warning('Protein not found for accession {}'.format(accession))
                        continue

                    if p:
                        # create a link
                        link, created = WebLink.objects.get_or_create(web_resource=resource, index=link_index)
                        if created:
                            num_created_links += 1

                        # add the link to this protein
                        p.web_links.add(link)

            if num_created_links:
                self.logger.info('Created {} web links for resource {}'.format(num_created_links, resource.name))
        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('COMPLETED CREATING LINKS')
