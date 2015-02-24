# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='webresource',
            name='slug',
            field=models.SlugField(default='', max_length=20),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='webresource',
            name='name',
            field=models.CharField(default='', max_length=200),
            preserve_default=True,
        ),
    ]
