from django.db import models
from string import Template


class WebResource(models.Model):
    slug = models.CharField(max_length=20)
    name = models.CharField(max_length=200, default='')
    url = models.TextField()
    #url should be a string template, so it can be automaticaly filled with index in proper place
    #https://docs.python.org/3.4/library/string.html#string.Template
    #Example: 'http://www.ncbi.nlm.nih.gov/pubmed/$index'

    def __str__(self):
        return self.name
    
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