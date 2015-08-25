# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0001_initial'),
        ('ligand', '0002_auto_20150429_1202'),
    ]

    operations = [
        migrations.CreateModel(
            name='LigandProperities',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('smiles', models.TextField(null=True)),
                ('inchikey', models.CharField(max_length=50, unique=True, null=True)),
                ('ligand_type', models.ForeignKey(to='ligand.LigandType', null=True)),
                ('web_links', models.ManyToManyField(to='common.WebLink')),
            ],
            options={
                'db_table': 'ligand_properities',
            },
        ),
        migrations.RemoveField(
            model_name='ligandalias',
            name='ligand',
        ),
        migrations.RemoveField(
            model_name='ligand',
            name='inchikey',
        ),
        migrations.RemoveField(
            model_name='ligand',
            name='ligand_type',
        ),
        migrations.RemoveField(
            model_name='ligand',
            name='smiles',
        ),
        migrations.RemoveField(
            model_name='ligand',
            name='web_links',
        ),
        migrations.AddField(
            model_name='ligand',
            name='ambigious_alias',
            field=models.NullBooleanField(),
        ),
        migrations.AddField(
            model_name='ligand',
            name='canonical',
            field=models.NullBooleanField(),
        ),
        migrations.DeleteModel(
            name='LigandAlias',
        ),
        migrations.AddField(
            model_name='ligand',
            name='properities',
            field=models.ForeignKey(default='', to='ligand.LigandProperities'),
            preserve_default=False,
        ),
    ]
