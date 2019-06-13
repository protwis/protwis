from django.db import models
from django.db import IntegrityError
from django.utils.text import slugify

from common.tools import fetch_from_web_api, fetch_from_entrez

from Bio import Entrez, Medline
from string import Template
import urllib.request,json
import logging
import re


class WebResource(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=200, default='')
    url = models.TextField()
    # url should be a string template, so it can be automaticaly filled with index in proper place
    # https://docs.python.org/3.4/library/string.html#string.Template
    # Example: 'http://www.ncbi.nlm.nih.gov/pubmed/$index'

    def __str__(self):
        return self.url
    
    class Meta():
        db_table = 'web_resource'

class WebLink(models.Model):
    web_resource = models.ForeignKey('WebResource', on_delete=models.CASCADE)
    index = models.TextField()

    # And now generating the working url is just a piece of cake
    def __str__(self):
        return Template(str(self.web_resource)).substitute(index=self.index)
    
    class Meta():
        db_table = 'web_link'
        unique_together = ('web_resource', 'index')


class Publication(models.Model):
    web_link = models.OneToOneField('common.WebLink', on_delete=models.CASCADE)
    journal = models.ForeignKey('PublicationJournal', null=True, on_delete=models.CASCADE)
    title = models.TextField(null=True)
    authors = models.TextField(null=True)
    year = models.IntegerField(null=True)
    reference = models.TextField(null=True)

    def __str__(self):
        return "{!s} ({!s})".format(self.journal, self.year)

    class Meta():
        db_table = 'publication'


    #http://www.ncbi.nlm.nih.gov/pubmed/?term=10.1124%2Fmol.107.040097&report=xml&format=text
    # use NCBI instead to correct year published (journal year)

    def update_from_doi(self, doi):
        logger = logging.getLogger('build')
        # should entrez be tried as a backup?
        try_entrez_on_fail = False
        
        # check whether this data is cached
        cache_dir = ['crossref', 'doi']
        url = 'http://api.crossref.org/works/$index'
        pub = fetch_from_web_api(url, doi, cache_dir)
                
        if pub:
            # update record
            try:
                self.title = pub['message']['title'][0]
                try: 
                    self.year = pub['message']['created']['date-parts'][0][0]
                except:
                    self.year = pub['message']['deposited']['date-parts'][0][0]

                # go from [{'family': 'Gloriam', 'given': 'David E.'}] to ['Gloriam DE']
                authors = ['{} {}'.format(x['family'], ''.join([y[:1] for y in x['given'].split()]))
                    for x in pub['message']['author']]
                self.authors = ', '.join(authors)
            
                # get volume and pages if available
                reference = {}
                fields = ['volume', 'page']
                for f in fields:
                    if f in pub['message']:
                        reference[f] = pub['message'][f]
                    else:
                        reference[f] = 'X'
                self.reference = '{}:{}'.format(reference['volume'], reference['page'])

                # journal
                journal = pub['message']['container-title'][0]
                try:
                    # not all records have the journal abbreviation
                    journal_abbr = pub['message']['container-title'][1]
                except:
                    journal_abbr = slugify(journal)
                try:
                    self.journal, created = PublicationJournal.objects.get_or_create(name=journal,
                        defaults={'slug': journal_abbr})
                    if created:
                        logger.info('Created journal {}'.format(journal))
                except IntegrityError:
                    self.journal = PublicationJournal.objects.get(name=journal)
            except Exception as msg:
                logger.warning('Processing data from CrossRef for {} failed: {}'.format(doi, msg))
                try_entrez_on_fail = False
        else:
            print("Publication not on crossref",doi)
            try_entrez_on_fail = False

        if try_entrez_on_fail:
            # try searching entrez for DOI
            try:
                Entrez.email = 'info@gpcrdb.org'
                record = Entrez.read(Entrez.esearch(
                    db='pubmed',
                    retmax=1,
                    term=doi
                    ))
                self.update_from_pubmed_data(record['IdList'][0])
            except:
                return False

    def update_from_pubmed_data(self, index=None):
        logger = logging.getLogger('build')

        if not index:
            index = self.web_link.index
        cache_dir = ['entrez', 'pmid']
        record = fetch_from_entrez(index, cache_dir)
        try:
            self.title = record['TI']
            self.authors = ', '.join(record['AU'])
            try:
                self.year = record['DA'][:4]
            except:
                # Sometimes 'DA' field does not exist, use alternative
                self.year = record['DP'][:4]

            record['JT'] = record['JT']
            record['TA'] = record['TA']
            try:
                self.journal, created = PublicationJournal.objects.get_or_create(name=record['JT'],
                        defaults={'slug': record['TA']})
            except PublicationJournal.DoesNotExist:
                j = PublicationJournal(slug=record['TA'], name=record['JT'])
                j.save()
                self.journal = j
            
            self.reference = ""
            if 'VI' in record:
                self.reference += record['VI']
            if 'PG' in record:
                self.reference += ":" + record['PG']
        except Exception as msg:
            logger.warning('Publication update on pubmed error! Pubmed: {} error {}'.format(index, msg))


class PublicationJournal(models.Model):
    slug = models.CharField(max_length=200, null=True)
    name = models.TextField(unique=True)

    def __str__(self):
        return "{!s} ({!s})".format(self.name, self.slug)

    class Meta():
        db_table = 'publication_journal'


class ReleaseNotes(models.Model):
    date = models.DateField()
    html = models.TextField()

    def __str__(self):
        return str(self.date)

    class Meta():
        ordering = ('-date', )
        db_table = 'release_notes'


class ReleaseStatistics(models.Model):
    release = models.ForeignKey('ReleaseNotes', on_delete=models.CASCADE)
    statistics_type = models.ForeignKey('ReleaseStatisticsType', on_delete=models.CASCADE)
    value = models.IntegerField()

    def __str__(self):
        return "{} {}".format(self.date, self.statistics_type)

    class Meta():
        ordering = ('id', )
        db_table = 'release_statistics'


class ReleaseStatisticsType(models.Model):
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'release_statistics_type'

