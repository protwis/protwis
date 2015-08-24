# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0005_remove_residue_sequence_based_generic_number'),
    ]

    operations = [
        migrations.CreateModel(
            name='ResidueGenericNumberEquivalent',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', primary_key=True, serialize=False)),
                ('label', models.CharField(max_length=10, db_index=True)),
                ('default_generic_number', models.ForeignKey(to='residue.ResidueGenericNumber')),
                ('scheme', models.ForeignKey(to='residue.ResidueNumberingScheme')),
            ],
            options={
                'db_table': 'residue_generic_number_equivalent',
            },
        ),
    ]
