from django.core.management.base import BaseCommand, CommandError

from common.models import Publication

import logging
import re


class Command(BaseCommand):
    help = 'Add (or refresh) Publication entries from DOIs, e.g. for the reference buttons (refLink) of the menus. ' + \
    'Usage: python manage.py add_publication https://doi.org/10.1093/nar/gkag1185 [more DOIs] [--year 2027] [--force]'
    logger = logging.getLogger('build')

    def add_arguments(self, parser):
        parser.add_argument('doi', nargs='+',
            help='DOI(s), either bare (10.1093/nar/gkag1185) or as link (https://doi.org/10.1093/nar/gkag1185)')
        parser.add_argument('--year', type=int, default=None,
            help='Override the publication year (e.g. journal issue year instead of online date)')
        parser.add_argument('-f', '--force',
            action='store_true',
            dest='force',
            help='Refresh metadata of publications that already exist')

    def handle(self, *args, **options):
        failed = []
        for doi_input in options['doi']:
            doi = self.parse_DOI(doi_input)
            if not doi:
                failed.append(doi_input)
                continue

            pubs = Publication.objects.filter(web_link__index__iexact=doi, web_link__web_resource__slug='doi')
            if pubs.exists():
                pub = pubs.first()
                if options['force']:
                    pub.update_from_doi(pub.web_link.index)
                    self.stdout.write('Refreshed existing publication {}'.format(doi))
                else:
                    self.stdout.write('Publication {} already exists (use --force to refresh)'.format(doi))
            else:
                pub = Publication.get_or_create_from_doi(doi)
                self.stdout.write('Created publication {}'.format(doi))

            if options['year']:
                pub.year = options['year']
            pub.save()

            if not pub.title:
                self.stderr.write('WARNING: no metadata found for {} (not on PubMed/Crossref?)'.format(doi))
                failed.append(doi_input)
            self.stdout.write('  Index:   {}\n  Title:   {}\n  Authors: {}\n  Journal: {}\n  Year:    {}'.format(
                pub.web_link.index, pub.title, pub.authors, pub.journal.name if pub.journal else None, pub.year))

        if failed:
            raise CommandError('Failed to add complete publication data for: {}'.format(', '.join(failed)))

    def parse_DOI(self, doi):
        # Strip doi.org link prefixes and store upper case, matching build_citations and getReference()
        doi = re.sub(r'^(https?://)?(dx\.)?doi\.org/', '', doi.strip(), flags=re.IGNORECASE)
        doi = re.sub(r'^doi:\s*', '', doi, flags=re.IGNORECASE)
        if not doi.startswith('10.'):
            self.stderr.write('Not a valid DOI: {}'.format(doi))
            return None
        return doi.upper()
