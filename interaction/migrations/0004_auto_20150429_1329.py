# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0002_auto_20150429_1313'),
        ('interaction', '0003_residuefragmentatom'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='rotamer',
            field=models.ForeignKey(default=0, to='structure.Rotamer'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='sfragment',
            field=models.ForeignKey(default=0, to='structure.Fragment'),
            preserve_default=False,
        ),
    ]
