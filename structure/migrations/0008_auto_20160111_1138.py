# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0007_auto_20160111_1106'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structurestabilizingagent',
            name='slug',
            field=models.SlugField(unique=True),
        ),
    ]
