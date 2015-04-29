# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0005_residuefragmentatom_interaction'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuefragmentatom',
            name='atomclass',
            field=models.CharField(default='', max_length=20),
            preserve_default=False,
        ),
    ]
