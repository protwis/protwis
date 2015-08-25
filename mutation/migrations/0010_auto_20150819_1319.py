# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0003_auto_20150814_1448'),
        ('mutation', '0009_auto_20150811_1529'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='mutationexperiment',
            name='ligand_class',
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='ligand_role',
            field=models.ForeignKey(to='ligand.LigandRole', null=True),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand', null=True),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='ligand_ref',
            field=models.ForeignKey(to='mutation.MutationLigandRef', null=True),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='refs',
            field=models.ForeignKey(to='common.Publication', null=True),
        ),
    ]
