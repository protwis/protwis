# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0003_auto_20150414_1459'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ligand',
            name='ligand_type',
            field=models.ForeignKey(to='ligand.LigandType', null=True),
        ),
    ]
