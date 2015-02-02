# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelTable(
            name='residue',
            table='residue',
        ),
        migrations.AlterModelTable(
            name='residuenumber',
            table='residue_number',
        ),
        migrations.AlterModelTable(
            name='residueset',
            table='residue_set',
        ),
    ]
