# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0007_auto_20150416_1121'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='residue',
            options={'ordering': ['sequence_number']},
        ),
    ]
