# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0010_auto_20150710_1457'),
        ('ligand', '0002_auto_20150429_1202'),
    ]

    operations = [
        migrations.CreateModel(
            name='AuxProtein',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('name', models.CharField(blank=True, max_length=100)),
                ('uniprot_id', models.CharField(max_length=20)),
                ('sequence', models.TextField(blank=True)),
                ('deletions', models.TextField(blank=True, max_length=100)),
                ('position', models.TextField(blank=True, max_length=20)),
                ('remarks', models.TextField(null=True)),
            ],
            options={
                'db_table': 'aux_protein',
            },
        ),
        migrations.CreateModel(
            name='AuxProteinType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('name', models.TextField(max_length=50)),
            ],
            options={
                'db_table': 'aux_protein_type',
            },
        ),
        migrations.CreateModel(
            name='Chemical',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=200)),
            ],
            options={
                'db_table': 'chemical',
            },
        ),
        migrations.CreateModel(
            name='ChemicalConc',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('concentration', models.TextField(null=True)),
                ('chemical', models.ForeignKey(to='construct.Chemical')),
            ],
            options={
                'db_table': 'chemical_conc',
            },
        ),
        migrations.CreateModel(
            name='ChemicalList',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('chemicals', models.ManyToManyField(to='construct.ChemicalConc')),
            ],
            options={
                'db_table': 'chemical_list',
            },
        ),
        migrations.CreateModel(
            name='ChemicalModification',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('description', models.TextField()),
            ],
            options={
                'db_table': 'chemical_modification',
            },
        ),
        migrations.CreateModel(
            name='ChemicalType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'chemical_type',
            },
        ),
        migrations.CreateModel(
            name='Construct',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('deletions', models.TextField(max_length=100)),
                ('parent', models.ForeignKey(to='construct.Construct', null=True)),
                ('protein_conformation', models.ForeignKey(to='protein.ProteinConformation')),
            ],
            options={
                'db_table': 'construct',
            },
        ),
        migrations.CreateModel(
            name='ConstructCrystallization',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('method', models.TextField(max_length=100, null=True)),
                ('settings', models.TextField(max_length=100, null=True)),
                ('remarks', models.TextField(null=True)),
                ('protein_conc', models.SlugField(max_length=20, blank=True)),
                ('aqueous_solution_lipid_ratio', models.SlugField(max_length=20, null=True)),
                ('lcp_bolus_volume', models.SlugField(max_length=20, null=True)),
                ('precipitant_solution_volume', models.SlugField(max_length=20, null=True)),
                ('temp', models.CharField(max_length=5, null=True)),
                ('ph', models.TextField(max_length=10, null=True)),
                ('chemical_list', models.ForeignKey(to='construct.ChemicalList')),
                ('construct', models.ForeignKey(to='construct.Construct')),
            ],
            options={
                'db_table': 'construct_crystallization',
            },
        ),
        migrations.CreateModel(
            name='ConstructCrystallizationLigandConc',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('ligand_conc', models.TextField(null=True)),
                ('construct_crystallization', models.ForeignKey(to='construct.ConstructCrystallization')),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
            ],
            options={
                'db_table': 'ligand_conc_of_crystallization',
            },
        ),
        migrations.CreateModel(
            name='ConstructExpression',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('remarks', models.TextField(null=True)),
                ('construct', models.ForeignKey(to='construct.Construct')),
            ],
            options={
                'db_table': 'construct_expression',
            },
        ),
        migrations.CreateModel(
            name='ConstructExpressionSystem',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('expression_method', models.CharField(max_length=100)),
                ('host_cell_type', models.CharField(max_length=100)),
                ('host_cell', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'construct_expression_system',
            },
        ),
        migrations.CreateModel(
            name='ConstructPurification',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('remarks', models.TextField(null=True)),
                ('construct', models.ForeignKey(to='construct.Construct')),
            ],
            options={
                'db_table': 'construct_purification',
            },
        ),
        migrations.CreateModel(
            name='ConstructSolubilization',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('remarks', models.TextField(null=True)),
                ('chemical_list', models.ForeignKey(to='construct.ChemicalList')),
                ('construct', models.ForeignKey(to='construct.Construct')),
            ],
            options={
                'db_table': 'construct_solubilization',
            },
        ),
        migrations.CreateModel(
            name='CrystallizationMethodTypes',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=100, null=True)),
            ],
            options={
                'db_table': 'crystallization_method_types',
            },
        ),
        migrations.CreateModel(
            name='PurificationStep',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('description', models.TextField(null=True)),
                ('purification', models.ForeignKey(to='construct.ConstructPurification')),
            ],
            options={
                'db_table': 'purification_step',
            },
        ),
        migrations.CreateModel(
            name='PurificationStepType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('name', models.TextField()),
            ],
            options={
                'db_table': 'purification_step_type',
            },
        ),
        migrations.AddField(
            model_name='purificationstep',
            name='purification_type',
            field=models.ForeignKey(to='construct.PurificationStepType'),
        ),
        migrations.AddField(
            model_name='constructexpression',
            name='expression_system',
            field=models.ForeignKey(to='construct.ConstructExpressionSystem'),
        ),
        migrations.AddField(
            model_name='constructcrystallization',
            name='crystal_type',
            field=models.ForeignKey(to='construct.CrystallizationMethodTypes'),
        ),
        migrations.AddField(
            model_name='constructcrystallization',
            name='ligands',
            field=models.ManyToManyField(to='ligand.Ligand', through='construct.ConstructCrystallizationLigandConc'),
        ),
        migrations.AddField(
            model_name='chemicalmodification',
            name='construct_solubilization',
            field=models.ForeignKey(to='construct.ConstructSolubilization'),
        ),
        migrations.AddField(
            model_name='chemical',
            name='chemical_type',
            field=models.ForeignKey(to='construct.ChemicalType'),
        ),
        migrations.AddField(
            model_name='auxprotein',
            name='construct',
            field=models.ForeignKey(to='construct.Construct'),
        ),
        migrations.AddField(
            model_name='auxprotein',
            name='protein_type',
            field=models.ForeignKey(to='construct.AuxProteinType'),
        ),
    ]
