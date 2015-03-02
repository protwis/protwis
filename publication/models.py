from django.db import models
from Bio import Entrez, Medline

class Publication(models.Model):
    web_link = models.ForeignKey('common.WebLink')
    title = models.TextField()
    year = models.IntegerField()
    journal = models.ForeignKey('PublicationJournal')
    citation = models.TextField()

    def update_from_pubmed_data(self, index=None):
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

