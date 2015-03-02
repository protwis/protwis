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
            name='role',
            field=models.ForeignKey(null=True, to='ligand.LigandRole'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='ligand',
            name='smiles',
            field=models.TextField(null=True),
            preserve_default=True,
        ),
    ]
