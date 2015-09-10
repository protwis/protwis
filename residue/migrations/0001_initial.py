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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('sequence_number', models.SmallIntegerField()),
                ('amino_acid', models.CharField(max_length=1)),
            ],
            options={
                'db_table': 'residue',
                'ordering': ['sequence_number'],
            },
        ),
        migrations.CreateModel(
            name='ResidueGenericNumber',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('label', models.CharField(db_index=True, max_length=10)),
                ('protein_segment', models.ForeignKey(null=True, to='protein.ProteinSegment', related_name='generic_numbers')),
            ],
            options={
                'db_table': 'residue_generic_number',
            },
        ),
        migrations.CreateModel(
            name='ResidueGenericNumberEquivalent',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('label', models.CharField(db_index=True, max_length=10)),
                ('default_generic_number', models.ForeignKey(to='residue.ResidueGenericNumber')),
            ],
            options={
                'db_table': 'residue_generic_number_equivalent',
            },
        ),
        migrations.CreateModel(
            name='ResidueNumberingScheme',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('slug', models.SlugField(max_length=20)),
                ('short_name', models.CharField(max_length=20)),
                ('name', models.CharField(max_length=100)),
                ('parent', models.ForeignKey(null=True, to='residue.ResidueNumberingScheme')),
            ],
            options={
                'db_table': 'residue_generic_numbering_scheme',
            },
        ),
        migrations.CreateModel(
            name='ResidueSet',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('name', models.CharField(max_length=50)),
                ('residue', models.ManyToManyField(to='residue.Residue')),
            ],
            options={
                'db_table': 'residue_set',
            },
        ),
        migrations.AddField(
            model_name='residuegenericnumberequivalent',
            name='scheme',
            field=models.ForeignKey(to='residue.ResidueNumberingScheme'),
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
            field=models.ForeignKey(null=True, to='residue.ResidueGenericNumber', related_name='display'),
        ),
        migrations.AddField(
            model_name='residue',
            name='generic_number',
            field=models.ForeignKey(null=True, to='residue.ResidueGenericNumber', related_name='compare'),
        ),
        migrations.AddField(
            model_name='residue',
            name='protein_conformation',
            field=models.ForeignKey(to='protein.ProteinConformation'),
        ),
        migrations.AddField(
            model_name='residue',
            name='protein_segment',
            field=models.ForeignKey(null=True, to='protein.ProteinSegment'),
        ),
    ]
