# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='structure',
            name='structure_type',
            field=models.CharField(null=True, max_length=20),
            preserve_default=True,
        ),
    ]
