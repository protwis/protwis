# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0007_protein_protein_anomalies'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='proteinanomaly',
            options={'ordering': ('generic_number__label',)},
        ),
    ]
