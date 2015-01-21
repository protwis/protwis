# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0002_auto_20150121_1621'),
        ('ligand', '0001_initial'),
        ('protein', '0002_auto_20150121_1621'),
        ('crystalstructure', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Fragment',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('filename', models.CharField(max_length=50)),
                ('crystalstructure', models.ForeignKey(to='crystalstructure.CrystalStructure')),
            ],
            options={
                'db_table': 'fragment',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Interaction',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('interaction_type', models.CharField(max_length=10)),
                ('description', models.CharField(max_length=200)),
            ],
            options={
                'db_table': 'interaction',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='fragment',
            name='interaction',
            field=models.ForeignKey(to='fragment.Interaction'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='fragment',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='fragment',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='fragment',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
            preserve_default=True,
        ),
    ]
