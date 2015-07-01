# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Residue',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('sequence_number', models.SmallIntegerField()),
                ('amino_acid', models.CharField(max_length=1)),
            ],
            options={
                'ordering': ['sequence_number'],
                'db_table': 'residue',
            },
        ),
        migrations.CreateModel(
            name='ResidueGenericNumber',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('label', models.CharField(max_length=10, db_index=True)),
                ('protein_segment', models.ForeignKey(to='protein.ProteinSegment', null=True)),
            ],
            options={
                'db_table': 'residue_generic_number',
            },
        ),
        migrations.CreateModel(
            name='ResidueNumberingScheme',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('short_name', models.CharField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'residue_generic_numbering_scheme',
            },
        ),
        migrations.CreateModel(
            name='ResidueSet',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('residue', models.ManyToManyField(to='residue.Residue')),
            ],
            options={
                'db_table': 'residue_set',
            },
        ),
        migrations.AddField(
            model_name='residuegenericnumber',
            name='scheme',
            field=models.ForeignKey(to='residue.ResidueNumberingScheme'),
        ),
        migrations.AddField(
            model_name='residue',
            name='alternative_generic_numbers',
            field=models.ManyToManyField(to='residue.ResidueGenericNumber', related_name='alternative'),
        ),
        migrations.AddField(
            model_name='residue',
            name='display_generic_number',
            field=models.ForeignKey(to='residue.ResidueGenericNumber', null=True, related_name='display'),
        ),
        migrations.AddField(
            model_name='residue',
            name='generic_number',
            field=models.ForeignKey(to='residue.ResidueGenericNumber', null=True, related_name='compare'),
        ),
        migrations.AddField(
            model_name='residue',
            name='protein_conformation',
            field=models.ForeignKey(to='protein.ProteinConformation'),
        ),
        migrations.AddField(
            model_name='residue',
            name='protein_segment',
            field=models.ForeignKey(to='protein.ProteinSegment', null=True),
        ),
    ]
