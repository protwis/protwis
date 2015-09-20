from django.db import models

# Create your models here.

class Documentation(models.Model):
    #Basic set 
    # defines the type of data
    title = models.TextField()
    description = models.TextField()
    image = models.TextField()
    html = models.TextField()

    def __str__(self):
        return self.title

    class Meta():
        db_table = "documentation"