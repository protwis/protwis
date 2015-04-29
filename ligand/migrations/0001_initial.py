# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Ligand',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.TextField()),
                ('smiles', models.TextField(null=True)),
                ('inchikey', models.CharField(max_length=50, null=True)),
            ],
            options={
                'db_table': 'ligand',
            },
        ),
        migrations.CreateModel(
            name='LigandAlias',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.TextField()),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
            ],
            options={
                'db_table': 'ligand_alias',
            },
        ),
        migrations.CreateModel(
            name='LigandRole',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField()),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'ligand_role',
            },
        ),
        migrations.CreateModel(
            name='LigandType',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField()),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'ligand_type',
            },
        ),
        migrations.AddField(
            model_name='ligand',
            name='ligand_type',
            field=models.ForeignKey(to='ligand.LigandType', null=True),
        ),
        migrations.AddField(
            model_name='ligand',
            name='web_links',
            field=models.ManyToManyField(to='common.WebLink'),
        ),
    ]
