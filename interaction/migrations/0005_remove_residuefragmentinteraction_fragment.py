# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0004_auto_20150429_1329'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='residuefragmentinteraction',
            name='fragment',
        ),
    ]
