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
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('sequence_number', models.SmallIntegerField()),
                ('generic_number', models.CharField(blank=True, max_length=10)),
                ('generic_number_alt', models.CharField(blank=True, max_length=10)),
                ('generic_number_alt2', models.CharField(blank=True, max_length=10)),
                ('generic_number_alt3', models.CharField(blank=True, max_length=10)),
                ('generic_number_alt4', models.CharField(blank=True, max_length=10)),
                ('amino_acid', models.CharField(max_length=1)),
                ('protein', models.ForeignKey(to='protein.Protein')),
                ('protein_segment', models.ForeignKey(to='protein.ProteinSegment', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ResidueNumber',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('generic_number', models.CharField(max_length=10)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ResidueSet',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('residue', models.ManyToManyField(to='residue.Residue')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
