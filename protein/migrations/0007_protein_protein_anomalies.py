# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0006_merge'),
    ]

    operations = [
        migrations.AddField(
            model_name='protein',
            name='protein_anomalies',
            field=models.ManyToManyField(to='protein.ProteinAnomaly'),
        ),
    ]
