# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2017-09-20 07:35
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='structure',
            name='refined',
            field=models.BooleanField(default=False),
        ),
    ]
