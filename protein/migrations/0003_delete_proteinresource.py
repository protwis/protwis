# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0002_auto_20150219_1527'),
    ]

    operations = [
        migrations.DeleteModel(
            name='ProteinResource',
        ),
    ]
