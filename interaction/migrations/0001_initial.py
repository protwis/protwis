# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0005_auto_20150408_1335'),
        ('protein', '0003_auto_20150412_2031'),
        ('ligand', '0003_auto_20150414_1459'),
        ('structure', '0003_auto_20150414_1459'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProteinLigandInteraction',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'db_table': 'interaction_protein_ligand',
            },
        ),
        migrations.CreateModel(
            name='ResidueFragmentInteraction',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('fragment', models.TextField()),
            ],
            options={
                'db_table': 'interaction_residue_fragment',
            },
        ),
        migrations.CreateModel(
            name='ResidueFragmentInteractionType',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
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
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
                ('structure', models.ForeignKey(to='structure.Structure')),
            ],
            options={
                'db_table': 'interaction_structure_ligand',
            },
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='interaction_type',
            field=models.ForeignKey(to='interaction.ResidueFragmentInteractionType'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
    ]
