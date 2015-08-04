# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0004_auto_20150727_1406'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mutation',
            name='foldchange',
            field=models.FloatField(),
        ),
    ]
