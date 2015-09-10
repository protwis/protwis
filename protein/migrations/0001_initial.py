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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('name', models.CharField(max_length=100)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'gene',
                'ordering': ('position',),
            },
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('entry_name', models.SlugField(unique=True, max_length=100)),
                ('accession', models.CharField(null=True, db_index=True, max_length=100)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('name', models.CharField(max_length=200)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'protein_alias',
                'ordering': ('position',),
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomaly',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
            ],
            options={
                'db_table': 'protein_anomaly',
                'ordering': ('generic_number__label',),
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyRule',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('amino_acid', models.CharField(max_length=1)),
                ('negative', models.BooleanField(default=False)),
            ],
            options={
                'db_table': 'protein_anomaly_rule',
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyRuleSet',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('exclusive', models.BooleanField(default=False)),
            ],
            options={
                'db_table': 'protein_anomaly_rule_set',
                'ordering': ('id',),
            },
        ),
        migrations.CreateModel(
            name='ProteinAnomalyType',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
            ],
            options={
                'db_table': 'protein_conformation',
            },
        ),
        migrations.CreateModel(
            name='ProteinConformationTemplateStructure',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('protein_conformation', models.ForeignKey(to='protein.ProteinConformation')),
            ],
            options={
                'db_table': 'protein_conformation_template_structure',
            },
        ),
        migrations.CreateModel(
            name='ProteinFamily',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('slug', models.SlugField(max_length=100)),
                ('name', models.CharField(max_length=200)),
                ('parent', models.ForeignKey(null=True, to='protein.ProteinFamily')),
            ],
            options={
                'db_table': 'protein_family',
                'ordering': ('id',),
            },
        ),
        migrations.CreateModel(
            name='ProteinFusion',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('slug', models.SlugField(max_length=100)),
                ('name', models.CharField(max_length=50)),
                ('category', models.CharField(max_length=50)),
                ('partial', models.BooleanField(default=False)),
            ],
            options={
                'db_table': 'protein_segment',
                'ordering': ('id',),
            },
        ),
        migrations.CreateModel(
            name='ProteinSequenceType',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('name', models.CharField(max_length=50)),
                ('proteins', models.ManyToManyField(to='protein.Protein')),
            ],
            options={
                'db_table': 'protein_set',
            },
        ),
        migrations.CreateModel(
            name='ProteinSource',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'protein_source',
            },
        ),
        migrations.CreateModel(
            name='ProteinState',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('latin_name', models.CharField(max_length=100)),
                ('common_name', models.CharField(max_length=100, blank=True)),
            ],
            options={
                'db_table': 'species',
            },
        ),
        migrations.AddField(
            model_name='proteinfusionprotein',
            name='segment_after',
            field=models.ForeignKey(related_name='segment_after', to='protein.ProteinSegment'),
        ),
        migrations.AddField(
            model_name='proteinfusionprotein',
            name='segment_before',
            field=models.ForeignKey(related_name='segment_before', to='protein.ProteinSegment'),
        ),
        migrations.AddField(
            model_name='proteinfusion',
            name='proteins',
            field=models.ManyToManyField(to='protein.Protein', through='protein.ProteinFusionProtein'),
        ),
        migrations.AddField(
            model_name='proteinconformationtemplatestructure',
            name='protein_segment',
            field=models.ForeignKey(to='protein.ProteinSegment'),
        ),
    ]
