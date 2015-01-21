# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Ligand',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('name', models.CharField(max_length=50)),
                ('xray_name', models.CharField(null=True, max_length=3)),
                ('smiles', models.TextField()),
            ],
            options={
                'db_table': 'ligand',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='LigandRole',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('role', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'ligand_role',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='ligand',
            name='role',
            field=models.ForeignKey(to='ligand.LigandRole'),
            preserve_default=True,
        ),
    ]
