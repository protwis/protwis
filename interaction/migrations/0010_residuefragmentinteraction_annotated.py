# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0009_auto_20150505_1451'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='annotated',
            field=models.BooleanField(default=False),
        ),
    ]
