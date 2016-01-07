# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0005_auto_20150915_0922'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='proteinanomaly',
            unique_together=set([('anomaly_type', 'generic_number')]),
        ),
    ]
