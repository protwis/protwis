# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0002_residue_sequence_based_generic_number'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuenumberingscheme',
            name='parent',
            field=models.ForeignKey(to='residue.ResidueNumberingScheme', null=True),
        ),
    ]
