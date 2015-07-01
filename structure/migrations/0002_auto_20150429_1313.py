# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('ligand', '0002_auto_20150429_1202'),
        ('structure', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Fragment',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('ligand', models.ForeignKey(to='ligand.Ligand')),
            ],
            options={
                'db_table': 'structure_fragment',
            },
        ),
        migrations.CreateModel(
            name='PdbData',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('pdb', models.TextField()),
            ],
            options={
                'db_table': 'structure_pdb_data',
            },
        ),
        migrations.CreateModel(
            name='Rotamer',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('pdbdata', models.ForeignKey(to='structure.PdbData')),
                ('residue', models.ForeignKey(to='residue.Residue')),
                ('structure', models.ForeignKey(to='structure.Structure')),
            ],
            options={
                'db_table': 'structure_rotamer',
            },
        ),
        migrations.AddField(
            model_name='fragment',
            name='pdbdata',
            field=models.ForeignKey(to='structure.PdbData'),
        ),
        migrations.AddField(
            model_name='fragment',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
    ]
