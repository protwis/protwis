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
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
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
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=200)),
                ('position', models.SmallIntegerField()),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'ordering': ('position',),
                'db_table': 'protein_alias',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinFamily',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(unique=True, max_length=100)),
                ('name', models.CharField(max_length=200)),
                ('parent', models.ForeignKey(null=True, to='protein.ProteinFamily')),
            ],
            options={
                'db_table': 'protein_family',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinLinks',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('url', models.TextField()),
            ],
            options={
                'db_table': 'protein_links',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinResource',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=200)),
                ('url', models.TextField()),
            ],
            options={
                'db_table': 'protein_resource',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinSegment',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(unique=True, max_length=100)),
                ('name', models.CharField(max_length=50)),
                ('category', models.CharField(max_length=50)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'ordering': ('position',),
                'db_table': 'protein_segment',
            },
            bases=(models.Model,),
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
            bases=(models.Model,),
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
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('latin_name', models.CharField(max_length=100)),
                ('common_name', models.CharField(max_length=100, blank=True)),
            ],
            options={
                'db_table': 'species',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='proteinlinks',
            name='resource',
            field=models.ForeignKey(to='protein.ProteinResource'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='protein',
            name='family',
            field=models.ForeignKey(to='protein.ProteinFamily'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='protein',
            name='source',
            field=models.ForeignKey(to='protein.ProteinSource'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='protein',
            name='species',
            field=models.ForeignKey(to='protein.Species'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='gene',
            name='proteins',
            field=models.ManyToManyField(to='protein.Protein'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='gene',
            name='species',
            field=models.ForeignKey(to='protein.Species'),
            preserve_default=True,
        ),
    ]
