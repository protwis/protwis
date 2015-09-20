# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0002_structure_representative'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structuretype',
            name='slug',
            field=models.SlugField(max_length=20, unique=True),
        ),
    ]
