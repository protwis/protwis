# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0008_auto_20150602_1321'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='protein',
            name='protein_anomalies',
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='protein_anomalies',
            field=models.ManyToManyField(to='protein.ProteinAnomaly'),
        ),
    ]
