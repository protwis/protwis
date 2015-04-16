# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0006_aminoacid_aminoacidset'),
        ('protein', '0004_protein_endogenous_ligand'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProteinAnomaly',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'db_table': 'protein_anomaly',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyRule',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('amino_acid', models.ManyToManyField(to='residue.AminoAcid')),
                ('generic_number', models.ForeignKey(to='residue.ResidueGenericNumber')),
            ],
            options={
                'db_table': 'protein_anomaly_rule',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyRuleSet',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('protein_anomaly', models.ForeignKey(to='protein.ProteinAnomaly')),
                ('protein_anomaly_rules', models.ManyToManyField(to='protein.ProteinAnomalyRule')),
            ],
            options={
                'db_table': 'protein_anomaly_rule_set',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyType',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('slug', models.SlugField(max_length=20, unique=True)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'protein_anomaly_type',
            },
        ),
        migrations.AddField(
            model_name='proteinanomaly',
            name='anomaly_type',
            field=models.ForeignKey(to='protein.ProteinAnomalyType'),
        ),
        migrations.AddField(
            model_name='proteinanomaly',
            name='generic_number',
            field=models.ForeignKey(to='residue.ResidueGenericNumber'),
        ),
    ]
