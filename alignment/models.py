from django.db import models

# Create your models here.

class AlignmentConsensus(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    alignment = models.BinaryField()