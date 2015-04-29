# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0005_auto_20150417_1117'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ligand',
            name='inchikey',
            field=models.CharField(unique=True, max_length=50, null=True),
        ),
        migrations.AlterField(
            model_name='ligand',
            name='ligand_type',
            field=models.ForeignKey(to='ligand.LigandType'),
        ),
        migrations.AlterField(
            model_name='ligandrole',
            name='slug',
            field=models.SlugField(max_length=10),
        ),
        migrations.AlterField(
            model_name='ligandtype',
            name='slug',
            field=models.SlugField(max_length=20, unique=True),
        ),
    ]
