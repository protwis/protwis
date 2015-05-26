# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0003_residuenumberingscheme_parent'),
    ]

    operations = [
        migrations.AlterField(
            model_name='residuegenericnumber',
            name='protein_segment',
            field=models.ForeignKey(related_name='generic_numbers', null=True, to='protein.ProteinSegment'),
        ),
    ]
