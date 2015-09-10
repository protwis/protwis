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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
            ],
            options={
                'db_table': 'interaction_protein_ligand',
            },
        ),
        migrations.CreateModel(
            name='ResidueFragmentAtom',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('atomtype', models.CharField(max_length=20)),
                ('atomnr', models.SmallIntegerField()),
                ('atomclass', models.CharField(max_length=20)),
                ('residuename', models.CharField(max_length=20)),
                ('chain', models.CharField(max_length=20)),
                ('residuenr', models.SmallIntegerField()),
                ('x', models.DecimalField(max_digits=6, decimal_places=3)),
                ('y', models.DecimalField(max_digits=6, decimal_places=3)),
                ('z', models.DecimalField(max_digits=6, decimal_places=3)),
                ('occupancy', models.DecimalField(max_digits=6, decimal_places=2)),
                ('temperature', models.DecimalField(max_digits=6, decimal_places=2)),
                ('element_name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'interaction_residue_fragment_atoms',
            },
        ),
        migrations.CreateModel(
            name='ResidueFragmentInteraction',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
            ],
            options={
                'db_table': 'interaction_residue_fragment',
            },
        ),
        migrations.CreateModel(
            name='ResidueFragmentInteractionType',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('pdb_reference', models.CharField(null=True, max_length=3)),
                ('annotated', models.BooleanField(default=False)),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
                ('ligand_role', models.ForeignKey(to='ligand.LigandRole')),
            ],
            options={
                'db_table': 'interaction_structure_ligand',
            },
        ),
    ]
