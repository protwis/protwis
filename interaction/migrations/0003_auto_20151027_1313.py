# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0002_auto_20150825_1114'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuefragmentinteractiontype',
            name='direction',
            field=models.CharField(null=True, max_length=30),
        ),
        migrations.AddField(
            model_name='residuefragmentinteractiontype',
            name='type',
            field=models.CharField(null=True, max_length=50),
        ),
    ]
