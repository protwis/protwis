# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0003_auto_20150414_1215'),
        ('ligand', '0002_auto_20150225_1441'),
    ]

    operations = [
        migrations.CreateModel(
            name='LigandAlias',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('name', models.TextField()),
            ],
            options={
                'db_table': 'ligand_alias',
            },
        ),
        migrations.CreateModel(
            name='LigandType',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'ligand_type',
            },
        ),
        migrations.RemoveField(
            model_name='ligand',
            name='role',
        ),
        migrations.RemoveField(
            model_name='ligand',
            name='xray_name',
        ),
        migrations.RemoveField(
            model_name='ligandrole',
            name='role',
        ),
        migrations.AddField(
            model_name='ligand',
            name='inchikey',
            field=models.CharField(null=True, max_length=50),
        ),
        migrations.AddField(
            model_name='ligand',
            name='web_link',
            field=models.ManyToManyField(to='common.WebLink'),
        ),
        migrations.AddField(
            model_name='ligandrole',
            name='name',
            field=models.CharField(default='', max_length=100),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='ligandrole',
            name='slug',
            field=models.SlugField(default='', max_length=10),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='ligand',
            name='name',
            field=models.TextField(),
        ),
        migrations.AddField(
            model_name='ligandalias',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand'),
        ),
        migrations.AddField(
            model_name='ligand',
            name='ligand_type',
            field=models.ForeignKey(default=0, to='ligand.LigandType'),
            preserve_default=False,
        ),
    ]
