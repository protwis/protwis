# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0001_initial'),
        ('protein', '0002_auto_20150219_1527'),
    ]

    operations = [
        migrations.CreateModel(
            name='Structure',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('pdb_code', models.CharField(max_length=4)),
                ('preferred_chain', models.CharField(max_length=20)),
                ('resolution', models.DecimalField(max_digits=5, decimal_places=3)),
                ('pmid', models.CharField(max_length=20)),
                ('publication_date', models.DateField()),
                ('n_terminus', models.TextField()),
                ('icl1', models.TextField()),
                ('ecl1', models.TextField()),
                ('icl2', models.TextField()),
                ('ecl2_1', models.TextField()),
                ('ecl2_2', models.TextField()),
                ('icl3', models.TextField()),
                ('ecl3', models.TextField()),
                ('c_terminus', models.TextField()),
                ('endogenous_ligand', models.ForeignKey(null=True, related_name='endogenous_ligand', to='ligand.Ligand')),
                ('protein_class', models.ForeignKey(to='protein.ProteinFamily')),
                ('receptor', models.ForeignKey(to='protein.Protein')),
                ('xray_ligand', models.ForeignKey(null=True, related_name='xray_ligand', to='ligand.Ligand')),
            ],
            options={
                'db_table': 'structure',
            },
            bases=(models.Model,),
        ),
    ]
