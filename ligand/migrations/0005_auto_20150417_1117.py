# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0004_auto_20150417_1042'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ligandrole',
            name='slug',
            field=models.SlugField(),
        ),
        migrations.AlterField(
            model_name='ligandtype',
            name='slug',
            field=models.SlugField(),
        ),
    ]
