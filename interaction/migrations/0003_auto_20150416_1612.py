# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0002_auto_20150414_2102'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structureligandinteraction',
            name='pdb_reference',
            field=models.CharField(null=True, max_length=3),
        ),
    ]
