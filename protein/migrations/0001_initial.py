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
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('name', models.CharField(max_length=100)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'ordering': ('position',),
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('accession', models.CharField(max_length=100)),
                ('entry_name', models.SlugField(max_length=100, unique=True)),
                ('name', models.CharField(max_length=200)),
                ('sequence', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinAlias',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('name', models.CharField(max_length=200)),
                ('position', models.SmallIntegerField()),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'ordering': ('position',),
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinFamily',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('slug', models.SlugField(max_length=100, unique=True)),
                ('name', models.CharField(max_length=200)),
                ('parent', models.ForeignKey(to='protein.ProteinFamily', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinLinks',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('url', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinResource',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('name', models.CharField(max_length=200)),
                ('url', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinSegment',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('slug', models.SlugField(max_length=100, unique=True)),
                ('name', models.CharField(max_length=50)),
                ('category', models.CharField(max_length=50)),
                ('position', models.SmallIntegerField()),
            ],
            options={
                'ordering': ('position',),
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinSet',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('protein', models.ManyToManyField(to='protein.Protein')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinSource',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('name', models.CharField(max_length=20)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('latin_name', models.CharField(max_length=100)),
                ('common_name', models.CharField(blank=True, max_length=100)),
            ],
            options={
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
