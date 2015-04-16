# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0005_auto_20150408_1335'),
    ]

    operations = [
        migrations.CreateModel(
            name='AminoAcid',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('one_letter_code', models.CharField(max_length=1)),
                ('three_letter_code', models.CharField(max_length=3)),
            ],
            options={
                'db_table': 'amino_acid',
            },
        ),
        migrations.CreateModel(
            name='AminoAcidSet',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('name', models.CharField(max_length=100)),
                ('amino_acids', models.ManyToManyField(to='residue.AminoAcid')),
            ],
            options={
                'db_table': 'amino_acid_set',
            },
        ),
    ]
