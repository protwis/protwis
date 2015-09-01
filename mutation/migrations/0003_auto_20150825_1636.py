# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0002_auto_20150825_1114'),
    ]

    operations = [
        migrations.DeleteModel(
            name='MutationLigand',
        ),
        migrations.DeleteModel(
            name='MutationRefs',
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='exp_func',
            field=models.ForeignKey(to='mutation.MutationFunc', null=True),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='exp_measure',
            field=models.ForeignKey(to='mutation.MutationMeasure', null=True),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='exp_qual',
            field=models.ForeignKey(to='mutation.MutationQual', null=True),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='exp_type',
            field=models.ForeignKey(to='mutation.MutationExperimentalType', null=True),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='optional',
            field=models.ForeignKey(to='mutation.MutationOptional', null=True),
        ),
    ]
