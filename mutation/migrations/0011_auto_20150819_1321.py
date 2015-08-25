# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0010_auto_20150819_1319'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mutationexperiment',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand', null=True, related_name='ligand'),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='ligand_ref',
            field=models.ForeignKey(to='ligand.Ligand', null=True, related_name='reference_ligand'),
        ),
    ]
