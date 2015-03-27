# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0003_auto_20150318_1550'),
    ]

    operations = [
        migrations.AlterField(
            model_name='residuegenericnumber',
            name='label',
            field=models.CharField(db_index=True, max_length=10),
            preserve_default=True,
        ),
    ]
