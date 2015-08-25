# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('common', '0001_initial'),
        ('mutation', '0001_initial'),
        ('protein', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='mutationexperiment',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='raw',
            field=models.ForeignKey(to='mutation.MutationRaw'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='refs',
            field=models.ForeignKey(null=True, to='common.Publication'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='mutation_type',
            field=models.ForeignKey(null=True, to='mutation.MutationType'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='residue',
            field=models.ForeignKey(null=True, to='residue.Residue'),
        ),
    ]
