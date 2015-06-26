# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0004_auto_20150515_1114'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='residue',
            name='sequence_based_generic_number',
        ),
    ]
