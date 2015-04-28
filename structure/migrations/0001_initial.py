# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0001_initial'),
        ('common', '0001_initial'),
        ('ligand', '0001_initial'),
        ('interaction', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Structure',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('preferred_chain', models.CharField(max_length=20)),
                ('resolution', models.DecimalField(decimal_places=3, max_digits=5)),
                ('publication_date', models.DateField()),
                ('ligands', models.ManyToManyField(to='ligand.Ligand', through='interaction.StructureLigandInteraction')),
                ('pdb_code', models.ForeignKey(to='common.WebLink')),
                ('protein_anomalies', models.ManyToManyField(to='protein.ProteinAnomaly')),
                ('protein_conformation', models.ForeignKey(to='protein.ProteinConformation')),
                ('publication', models.ForeignKey(to='common.Publication', null=True)),
            ],
            options={
                'db_table': 'structure',
            },
        ),
        migrations.CreateModel(
            name='StructureModel',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'db_table': 'structure_model',
            },
        ),
        migrations.CreateModel(
            name='StructureStabilizingAgent',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'structure_stabilizing_agent',
            },
        ),
        migrations.CreateModel(
            name='StructureType',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'structure_type',
            },
        ),
        migrations.AddField(
            model_name='structure',
            name='stabilizing_agents',
            field=models.ManyToManyField(to='structure.StructureStabilizingAgent'),
        ),
        migrations.AddField(
            model_name='structure',
            name='state',
            field=models.ForeignKey(to='protein.ProteinState'),
        ),
        migrations.AddField(
            model_name='structure',
            name='structure_type',
            field=models.ForeignKey(to='structure.StructureType'),
        ),
    ]
