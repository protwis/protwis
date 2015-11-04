# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0003_auto_20151019_0850'),
    ]

    operations = [
        migrations.AlterField(
            model_name='publicationjournal',
            name='name',
            field=models.TextField(unique=True),
        ),
    ]
