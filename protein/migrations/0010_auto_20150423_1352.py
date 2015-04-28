# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0009_auto_20150422_0910'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProteinConformation',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', auto_created=True, serialize=False)),
            ],
            options={
                'db_table': 'protein_conformation',
            },
        ),
        migrations.CreateModel(
            name='ProteinState',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', auto_created=True, serialize=False)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'protein_state',
            },
        ),
        migrations.RenameField(
            model_name='protein',
            old_name='endogenous_ligand',
            new_name='endogenous_ligands',
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='state',
            field=models.ForeignKey(to='protein.ProteinState'),
        ),
        migrations.AddField(
            model_name='protein',
            name='states',
            field=models.ManyToManyField(through='protein.ProteinConformation', to='protein.ProteinState'),
        ),
    ]
