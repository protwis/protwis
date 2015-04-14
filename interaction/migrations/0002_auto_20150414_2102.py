# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0003_auto_20150414_1459'),
        ('interaction', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='structureligandinteraction',
            name='ligand_role',
            field=models.ForeignKey(to='ligand.LigandRole', default=0),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='structureligandinteraction',
            name='pdb_reference',
            field=models.CharField(max_length=3, default=''),
            preserve_default=False,
        ),
    ]
