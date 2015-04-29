# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0008_auto_20150428_1308'),
    ]

    operations = [
        migrations.AlterField(
            model_name='residuefragmentatom',
            name='occupancy',
            field=models.DecimalField(max_digits=6, decimal_places=2),
        ),
        migrations.AlterField(
            model_name='residuefragmentatom',
            name='temperature',
            field=models.DecimalField(max_digits=6, decimal_places=2),
        ),
    ]
