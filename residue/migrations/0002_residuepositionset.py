# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='ResiduePositionSet',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('residue_position', models.ManyToManyField(to='residue.ResidueGenericNumberEquivalent')),
            ],
            options={
                'db_table': 'residue_position_set',
            },
        ),
    ]
