from django.db import models

# Create your models here.

class News(models.Model):
    #Basic set 
    # defines the type of data
    image = models.TextField()
    date = models.DateField()
    html = models.TextField()

    def __str__(self):
        return self.title

    class Meta():
        db_table = "news"