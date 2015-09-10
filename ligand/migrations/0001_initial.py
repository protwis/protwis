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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('name', models.TextField()),
                ('canonical', models.NullBooleanField()),
                ('ambigious_alias', models.NullBooleanField()),
            ],
            options={
                'db_table': 'ligand',
            },
        ),
        migrations.CreateModel(
            name='LigandProperities',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('smiles', models.TextField(null=True)),
                ('inchikey', models.CharField(null=True, max_length=50, unique=True)),
            ],
            options={
                'db_table': 'ligand_properities',
            },
        ),
        migrations.CreateModel(
            name='LigandRole',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('slug', models.SlugField(unique=True, max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'ligand_type',
            },
        ),
        migrations.AddField(
            model_name='ligandproperities',
            name='ligand_type',
            field=models.ForeignKey(null=True, to='ligand.LigandType'),
        ),
        migrations.AddField(
            model_name='ligandproperities',
            name='web_links',
            field=models.ManyToManyField(to='common.WebLink'),
        ),
        migrations.AddField(
            model_name='ligand',
            name='properities',
            field=models.ForeignKey(to='ligand.LigandProperities'),
        ),
    ]
