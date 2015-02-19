from django.db import models


class WebResource(models.Model):
    name = models.CharField(max_length=200)
    url = models.TextField()

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'web_resource'

class WebLink(models.Model):
    web_resource = models.ForeignKey('WebResource')
    index = models.TextField()

    def __str__(self):
        return self.url
    
    class Meta():
        db_table = 'web_link'