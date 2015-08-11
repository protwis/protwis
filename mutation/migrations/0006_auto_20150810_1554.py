# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0010_auto_20150710_1457'),
        ('residue', '0005_remove_residue_sequence_based_generic_number'),
        ('mutation', '0005_auto_20150727_1443'),
    ]

    operations = [
        migrations.CreateModel(
            name='MutationExperiment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('mutation_to', models.CharField(max_length=1)),
                ('wt_value', models.DecimalField(max_digits=10, decimal_places=2)),
                ('wt_unit', models.CharField(max_length=10)),
                ('mu_value', models.DecimalField(max_digits=10, decimal_places=2)),
                ('mu_sign', models.CharField(max_length=2)),
                ('foldchange', models.FloatField()),
                ('exp_func', models.ForeignKey(to='mutation.MutationFunc')),
                ('exp_measure', models.ForeignKey(to='mutation.MutationMeasure')),
                ('exp_qual', models.ForeignKey(to='mutation.MutationQual')),
            ],
            options={
                'db_table': 'mutation_experiment',
            },
        ),
        migrations.CreateModel(
            name='MutationExperimentalType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_experimental_type',
            },
        ),
        migrations.RenameField(
            model_name='mutation',
            old_name='exp_type',
            new_name='mutation_type',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='exp_func',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='exp_measure',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='exp_qual',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='foldchange',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='ligand',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='ligand_class',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='ligand_ref',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='mu_sign',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='mu_value',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='optional',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='raw',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='refs',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='wt_unit',
        ),
        migrations.RemoveField(
            model_name='mutation',
            name='wt_value',
        ),
        migrations.AlterField(
            model_name='mutation',
            name='residue',
            field=models.ForeignKey(null=True, to='residue.Residue'),
        ),
        migrations.AlterModelTable(
            name='mutation',
            table=None,
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='exp_type',
            field=models.ForeignKey(to='mutation.MutationExperimentalType'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='ligand',
            field=models.ForeignKey(null=True, to='mutation.MutationLigand'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='ligand_class',
            field=models.ForeignKey(to='mutation.MutationLigandClass'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='ligand_ref',
            field=models.ForeignKey(to='mutation.MutationLigandRef'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='optional',
            field=models.ForeignKey(to='mutation.MutationOptional'),
        ),
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
            field=models.ForeignKey(null=True, to='mutation.MutationRefs'),
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
        ),
    ]
