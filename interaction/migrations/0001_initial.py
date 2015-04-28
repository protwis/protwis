# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProteinLigandInteraction',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'interaction_protein_ligand',
            },
        ),
        migrations.CreateModel(
            name='ResidueFragmentInteraction',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('fragment', models.TextField()),
            ],
            options={
                'db_table': 'interaction_residue_fragment',
            },
        ),
        migrations.CreateModel(
            name='ResidueFragmentInteractionType',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=10)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'interaction_type_residue_fragment',
            },
        ),
        migrations.CreateModel(
            name='StructureLigandInteraction',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('pdb_reference', models.CharField(max_length=3, null=True)),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
                ('ligand_role', models.ForeignKey(to='ligand.LigandRole')),
            ],
            options={
                'db_table': 'interaction_structure_ligand',
            },
        ),
    ]
