# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ligand',
            name='inchikey',
            field=models.CharField(unique=True, null=True, max_length=50),
        ),
        migrations.AlterField(
            model_name='ligandtype',
            name='slug',
            field=models.SlugField(unique=True, max_length=20),
        ),
    ]
