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
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=100)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'gene',
                'ordering': ('position',),
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('accession', models.CharField(max_length=100)),
                ('entry_name', models.SlugField(unique=True, max_length=100)),
                ('name', models.CharField(max_length=200)),
                ('sequence', models.TextField()),
            ],
            options={
                'db_table': 'protein',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinAlias',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=200)),
                ('position', models.SmallIntegerField()),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'db_table': 'protein_alias',
                'ordering': ('position',),
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinFamily',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('slug', models.SlugField(unique=True, max_length=100)),
                ('name', models.CharField(max_length=200)),
                ('parent', models.ForeignKey(to='protein.ProteinFamily', null=True)),
            ],
            options={
                'db_table': 'protein_family',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinSegment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('slug', models.SlugField(unique=True, max_length=100)),
                ('name', models.CharField(max_length=50)),
                ('category', models.CharField(max_length=50)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'db_table': 'protein_segment',
                'ordering': ('position',),
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinSet',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('protein', models.ManyToManyField(to='protein.Protein')),
            ],
            options={
                'db_table': 'protein_set',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinSource',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'protein_source',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('latin_name', models.CharField(max_length=100)),
                ('common_name', models.CharField(max_length=100, blank=True)),
            ],
            options={
                'db_table': 'species',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='protein',
            name='family',
            field=models.ForeignKey(to='protein.ProteinFamily'),
            preserve_default=True,
        ),
    ]
