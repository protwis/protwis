# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0007_auto_20150916_1308'),
    ]

    operations = [
        migrations.AddField(
            model_name='proteinsegment',
            name='fully_aligned',
            field=models.BooleanField(default=False),
        ),
    ]
