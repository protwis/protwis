# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='residue',
            name='sequence_based_generic_number',
            field=models.ForeignKey(null=True, related_name='seq_compare', to='residue.ResidueGenericNumber'),
        ),
    ]
