# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Mutation',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('mutation_to', models.CharField(max_length=1)),
                ('wt_value', models.DecimalField(decimal_places=2, max_digits=10)),
                ('wt_unit', models.CharField(max_length=10)),
                ('mu_value', models.DecimalField(decimal_places=2, max_digits=10)),
                ('mu_sign', models.CharField(max_length=2)),
                ('foldchange', models.DecimalField(decimal_places=2, max_digits=10)),
            ],
            options={
                'db_table': 'mutation',
            },
        ),
        migrations.CreateModel(
            name='MutationFunc',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('func', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_func',
            },
        ),
        migrations.CreateModel(
            name='MutationLigand',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('idtype', models.CharField(max_length=100)),
                ('name', models.CharField(max_length=100)),
                ('idid', models.CharField(max_length=100)),
                ('longseq', models.TextField()),
            ],
            options={
                'db_table': 'mutation_ligands',
            },
        ),
        migrations.CreateModel(
            name='MutationLigandClass',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('classname', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_ligand_class',
            },
        ),
        migrations.CreateModel(
            name='MutationLigandRef',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('reference', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_ligand_reference',
            },
        ),
        migrations.CreateModel(
            name='MutationMeasure',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('measure', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_measure',
            },
        ),
        migrations.CreateModel(
            name='MutationOptional',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('type', models.CharField(max_length=100)),
                ('wt', models.DecimalField(decimal_places=2, max_digits=10)),
                ('mu', models.DecimalField(decimal_places=2, max_digits=10)),
                ('sign', models.CharField(max_length=2)),
                ('percentage', models.DecimalField(decimal_places=2, max_digits=10)),
                ('qual', models.CharField(max_length=100)),
                ('agonist', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_opt',
            },
        ),
        migrations.CreateModel(
            name='MutationQual',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('qual', models.CharField(max_length=100)),
                ('prop', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_qual',
            },
        ),
        migrations.CreateModel(
            name='MutationRaw',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('reference', models.CharField(max_length=100)),
                ('protein', models.CharField(max_length=100)),
                ('mutation_pos', models.SmallIntegerField()),
                ('mutation_from', models.CharField(max_length=1)),
                ('mutation_to', models.CharField(max_length=1)),
                ('ligand_name', models.CharField(max_length=100)),
                ('ligand_idtype', models.CharField(max_length=100)),
                ('ligand_id', models.CharField(max_length=100)),
                ('ligand_class', models.CharField(max_length=100)),
                ('exp_type', models.CharField(max_length=100)),
                ('exp_func', models.CharField(max_length=100)),
                ('exp_wt_value', models.DecimalField(decimal_places=2, max_digits=10)),
                ('exp_wt_unit', models.CharField(max_length=10)),
                ('exp_mu_effect_type', models.CharField(max_length=100)),
                ('exp_mu_effect_sign', models.CharField(max_length=2)),
                ('exp_mu_effect_value', models.DecimalField(decimal_places=2, max_digits=10)),
                ('exp_mu_effect_qual', models.CharField(max_length=100)),
                ('exp_mu_effect_ligand_prop', models.CharField(max_length=100)),
                ('exp_mu_ligand_ref', models.CharField(max_length=100)),
                ('opt_type', models.CharField(max_length=100)),
                ('opt_wt', models.DecimalField(decimal_places=2, max_digits=10)),
                ('opt_mu', models.DecimalField(decimal_places=2, max_digits=10)),
                ('opt_sign', models.CharField(max_length=5)),
                ('opt_percentage', models.DecimalField(decimal_places=2, max_digits=10)),
                ('opt_qual', models.CharField(max_length=100)),
                ('opt_agonist', models.CharField(max_length=100)),
                ('added_by', models.CharField(max_length=100)),
                ('added_date', models.DateField()),
            ],
            options={
                'db_table': 'mutation_raw',
            },
        ),
        migrations.CreateModel(
            name='MutationRefs',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('year', models.SmallIntegerField()),
                ('journal', models.CharField(max_length=100)),
                ('title', models.TextField()),
                ('citation', models.TextField()),
                ('link', models.URLField()),
                ('ref_type', models.CharField(max_length=100)),
                ('reference', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_refs',
            },
        ),
        migrations.CreateModel(
            name='MutationType',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('type', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'mutation_type',
            },
        ),
        migrations.AddField(
            model_name='mutation',
            name='exp_func',
            field=models.ForeignKey(to='mutation.MutationFunc'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='exp_measure',
            field=models.ForeignKey(to='mutation.MutationMeasure'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='exp_qual',
            field=models.ForeignKey(to='mutation.MutationQual'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='exp_type',
            field=models.ForeignKey(to='mutation.MutationType'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='ligand',
            field=models.ForeignKey(to='mutation.MutationLigand'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='ligand_class',
            field=models.ForeignKey(to='mutation.MutationLigandClass'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='ligand_ref',
            field=models.ForeignKey(to='mutation.MutationLigandRef'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='optional',
            field=models.ForeignKey(to='mutation.MutationOptional'),
        ),
    ]
