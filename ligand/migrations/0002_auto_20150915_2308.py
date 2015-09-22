# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ligandrole',
            name='slug',
            field=models.SlugField(unique=True),
        ),
    ]
