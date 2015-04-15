# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0003_auto_20150414_1459'),
        ('interaction', '0002_auto_20150414_2102'),
        ('structure', '0003_auto_20150414_1459'),
    ]

    operations = [
        migrations.AddField(
            model_name='structure',
            name='ligands',
            field=models.ManyToManyField(to='ligand.Ligand', through='interaction.StructureLigandInteraction'),
        ),
    ]
