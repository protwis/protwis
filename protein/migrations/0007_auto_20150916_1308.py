# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0006_auto_20150915_1112'),
    ]

    operations = [
        migrations.AlterField(
            model_name='proteinanomalytype',
            name='slug',
            field=models.SlugField(unique=True, max_length=20),
        ),
    ]
