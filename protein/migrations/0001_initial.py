# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=100)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'ordering': ('position',),
                'db_table': 'gene',
            },
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('entry_name', models.SlugField(max_length=100, unique=True)),
                ('accession', models.CharField(max_length=100, db_index=True, null=True)),
                ('name', models.CharField(max_length=200)),
                ('sequence', models.TextField()),
            ],
            options={
                'db_table': 'protein',
            },
        ),
        migrations.CreateModel(
            name='ProteinAlias',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=200)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'ordering': ('position',),
                'db_table': 'protein_alias',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomaly',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'protein_anomaly',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyRule',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('amino_acid', models.CharField(max_length=1)),
            ],
            options={
                'db_table': 'protein_anomaly_rule',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyRuleSet',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('protein_anomaly', models.ForeignKey(to='protein.ProteinAnomaly')),
            ],
            options={
                'db_table': 'protein_anomaly_rule_set',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyType',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'protein_anomaly_type',
            },
        ),
        migrations.CreateModel(
            name='ProteinConformation',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'db_table': 'protein_conformation',
            },
        ),
        migrations.CreateModel(
            name='ProteinFamily',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=100)),
                ('name', models.CharField(max_length=200)),
                ('parent', models.ForeignKey(to='protein.ProteinFamily', null=True)),
            ],
            options={
                'ordering': ('id',),
                'db_table': 'protein_family',
            },
        ),
        migrations.CreateModel(
            name='ProteinFusion',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=100, unique=True)),
                ('sequence', models.TextField(null=True)),
            ],
            options={
                'db_table': 'protein_fusion',
            },
        ),
        migrations.CreateModel(
            name='ProteinFusionProtein',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('protein', models.ForeignKey(to='protein.Protein')),
                ('protein_fusion', models.ForeignKey(to='protein.ProteinFusion')),
            ],
            options={
                'db_table': 'protein_fusion_protein',
            },
        ),
        migrations.CreateModel(
            name='ProteinSegment',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=100)),
                ('name', models.CharField(max_length=50)),
                ('category', models.CharField(max_length=50)),
            ],
            options={
                'ordering': ('id',),
                'db_table': 'protein_segment',
            },
        ),
        migrations.CreateModel(
            name='ProteinSequenceType',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'protein_sequence_type',
            },
        ),
        migrations.CreateModel(
            name='ProteinSet',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('protein', models.ManyToManyField(to='protein.Protein')),
            ],
            options={
                'db_table': 'protein_set',
            },
        ),
        migrations.CreateModel(
            name='ProteinSource',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'protein_source',
            },
        ),
        migrations.CreateModel(
            name='ProteinState',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'protein_state',
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('latin_name', models.CharField(max_length=100)),
                ('common_name', models.CharField(blank=True, max_length=100)),
            ],
            options={
                'db_table': 'species',
            },
        ),
        migrations.AddField(
            model_name='proteinfusionprotein',
            name='segment_after',
            field=models.ForeignKey(to='protein.ProteinSegment', related_name='segment_after'),
        ),
        migrations.AddField(
            model_name='proteinfusionprotein',
            name='segment_before',
            field=models.ForeignKey(to='protein.ProteinSegment', related_name='segment_before'),
        ),
        migrations.AddField(
            model_name='proteinfusion',
            name='proteins',
            field=models.ManyToManyField(to='protein.Protein', through='protein.ProteinFusionProtein'),
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='state',
            field=models.ForeignKey(to='protein.ProteinState'),
        ),
    ]
