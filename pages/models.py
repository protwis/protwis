from django.db import models

# Create your models here.

class Pages(models.Model):
    #Basic set 
    title = models.TextField()
    html = models.TextField()

    def __str__(self):
        return self.title

    class Meta():
        db_table = "pages"