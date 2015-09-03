# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0002_auto_20150825_1114'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='proteinconformation',
            options={'ordering': ('id',)},
        ),
    ]
