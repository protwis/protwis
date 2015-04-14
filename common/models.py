from django.db import models
from Bio import Entrez, Medline
from string import Template


class WebResource(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=200, default='')
    url = models.TextField()
    #url should be a string template, so it can be automaticaly filled with index in proper place
    #https://docs.python.org/3.4/library/string.html#string.Template
    #Example: 'http://www.ncbi.nlm.nih.gov/pubmed/$index'

    def __str__(self):
        return self.url
    
    class Meta():
        db_table = 'web_resource'

class WebLink(models.Model):
    web_resource = models.ForeignKey('WebResource')
    index = models.TextField()

    #And now generating the working url is just a piece of cake
    def __str__(self):
        return Template(str(self.web_resource)).substitute(index=self.index)
    
    class Meta():
        db_table = 'web_link'


class Publication(models.Model):
    web_link = models.ForeignKey('common.WebLink')
    title = models.TextField()
    year = models.IntegerField()
    journal = models.ForeignKey('PublicationJournal')
    citation = models.TextField()

    def __str__(self):
        return "{!s} ({!s})".format(self.journal, self.year)

    class Meta():
        db_table = 'publication'

    def update_from_pubmed_data(self, index=None):
        try:
            Entrez.email = 'A.N.Other@example.com'
            if index:
                handle = Entrez.efetch(
                    db="pubmed", 
                    id=index, 
                    rettype="medline", 
                    retmode="text"
                    )
            else:
                handle = Entrez.efetch(
                    db="pubmed", 
                    id=self.web_link.index,
                    rettype="medline",
                    retmode="text"
                    )
        except Exception as e:
            print("Failed to retrieve data for pubmed id: {}".format(index))
            return
        try:
            record = Medline.read(handle)
            self.title = record['TI']
            self.year = record['DA'][0:3]
            try:
                self.journal = PublicationJournal.objects.get(slug=record['TA'], name=record['JT'])
            except PublicationJournal.DoesNotExist:
                j = PublicationJournal(slug=record['TA'], name=record['JT'])
                j.save()
                self.journal = j
            self.citation = "{}({}):{}".format(record['VI'], record['IP'], record['PG'])
        except Exception as msg:
            print(msg)

class PublicationJournal(models.Model):
    slug = models.CharField(max_length=30, null=True)
    name = models.TextField()

    def __str__(self):
        return "{!s} ({!s})".format(self.name, self.slug)

    class Meta():
        db_table = 'publication_journal'