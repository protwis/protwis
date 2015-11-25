# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0003_auto_20151027_1313'),
    ]

    operations = [
        migrations.AlterField(
            model_name='residuefragmentinteractiontype',
            name='slug',
            field=models.SlugField(max_length=20),
        ),
    ]
