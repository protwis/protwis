# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0001_initial'),
        ('residue', '0001_initial'),
        ('ligand', '0001_initial'),
        ('common', '0001_initial'),
        ('protein', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Fragment',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
            ],
            options={
                'db_table': 'structure_fragment',
            },
        ),
        migrations.CreateModel(
            name='PdbData',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('pdb', models.TextField()),
            ],
            options={
                'db_table': 'structure_pdb_data',
            },
        ),
        migrations.CreateModel(
            name='Rotamer',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('pdbdata', models.ForeignKey(to='structure.PdbData')),
                ('residue', models.ForeignKey(to='residue.Residue')),
            ],
            options={
                'db_table': 'structure_rotamer',
            },
        ),
        migrations.CreateModel(
            name='Structure',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('preferred_chain', models.CharField(max_length=20)),
                ('resolution', models.DecimalField(max_digits=5, decimal_places=3)),
                ('publication_date', models.DateField()),
                ('ligands', models.ManyToManyField(to='ligand.Ligand', through='interaction.StructureLigandInteraction')),
                ('pdb_code', models.ForeignKey(to='common.WebLink')),
                ('pdb_data', models.ForeignKey(null=True, to='structure.PdbData')),
                ('protein_anomalies', models.ManyToManyField(to='protein.ProteinAnomaly')),
                ('protein_conformation', models.ForeignKey(to='protein.ProteinConformation')),
                ('publication', models.ForeignKey(null=True, to='common.Publication')),
            ],
            options={
                'db_table': 'structure',
            },
        ),
        migrations.CreateModel(
            name='StructureModel',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'db_table': 'structure_model',
            },
        ),
        migrations.CreateModel(
            name='StructureSegment',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('start', models.IntegerField()),
                ('end', models.IntegerField()),
                ('protein_segment', models.ForeignKey(to='protein.ProteinSegment')),
                ('structure', models.ForeignKey(to='structure.Structure')),
            ],
            options={
                'db_table': 'structure_segment',
            },
        ),
        migrations.CreateModel(
            name='StructureStabilizingAgent',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
        migrations.AddField(
            model_name='rotamer',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='fragment',
            name='pdbdata',
            field=models.ForeignKey(to='structure.PdbData'),
        ),
        migrations.AddField(
            model_name='fragment',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
        ),
        migrations.AddField(
            model_name='fragment',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
    ]
