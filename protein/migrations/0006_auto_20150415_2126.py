# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0005_auto_20150415_2125'),
    ]

    operations = [
        migrations.AlterModelTable(
            name='proteinanomalyruleset',
            table='protein_anomaly_rule_set',
        ),
    ]
