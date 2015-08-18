# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0010_residuefragmentinteraction_annotated'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='residuefragmentinteraction',
            name='annotated',
        ),
        migrations.AddField(
            model_name='structureligandinteraction',
            name='annotated',
            field=models.BooleanField(default=False),
        ),
    ]
