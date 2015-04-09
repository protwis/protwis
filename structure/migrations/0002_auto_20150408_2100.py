# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structure',
            name='publication',
            field=models.ForeignKey(to='common.Publication', null=True),
            preserve_default=True,
        ),
    ]
