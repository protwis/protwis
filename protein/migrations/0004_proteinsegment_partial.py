# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0003_auto_20150505_1302'),
    ]

    operations = [
        migrations.AddField(
            model_name='proteinsegment',
            name='partial',
            field=models.BooleanField(default=False),
        ),
    ]
