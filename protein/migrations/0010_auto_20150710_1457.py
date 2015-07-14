# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0009_auto_20150605_1015'),
    ]

    operations = [
        migrations.RenameField(
            model_name='proteinset',
            old_name='protein',
            new_name='proteins',
        ),
    ]
