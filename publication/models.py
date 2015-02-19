from django.db import models

class Publication(models.Model):
    web_link = models.ForeignKey('common.WebLink')
    title = models.TextField()
    year = models.IntegerField()
    journal = models.ForeignKey('PublicatioJournal')
    citation = models.TextField()


class PublicationJournal(models.Model):
    slug = models.CharField(max_length=30)
    name = models.TextField()

    def __str__(self):
        return "{!s} ({!s})".format(self.name, self.slug)

    class Meta():
        db_table = 'publication_journal'

