# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0006_auto_20150416_1121'),
        ('residue', '0006_aminoacid_aminoacidset'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='aminoacidset',
            name='amino_acids',
        ),
        migrations.DeleteModel(
            name='AminoAcid',
        ),
        migrations.DeleteModel(
            name='AminoAcidSet',
        ),
    ]
