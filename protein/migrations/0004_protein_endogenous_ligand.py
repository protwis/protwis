# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0003_auto_20150414_1459'),
        ('protein', '0003_auto_20150412_2031'),
    ]

    operations = [
        migrations.AddField(
            model_name='protein',
            name='endogenous_ligand',
            field=models.ManyToManyField(to='ligand.Ligand'),
        ),
    ]
