# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0005_auto_20150415_2149'),
        ('structure', '0004_structure_ligands'),
    ]

    operations = [
        migrations.AddField(
            model_name='structure',
            name='protein_anomaly',
            field=models.ManyToManyField(to='protein.ProteinAnomaly'),
        ),
    ]
