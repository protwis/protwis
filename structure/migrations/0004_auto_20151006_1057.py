# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0003_auto_20150915_1738'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structurestabilizingagent',
            name='slug',
            field=models.SlugField(max_length=20, unique=True),
        ),
    ]
