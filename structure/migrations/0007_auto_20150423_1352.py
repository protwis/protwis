# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0006_auto_20150416_1534'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structure',
            name='state',
            field=models.ForeignKey(to='protein.ProteinState'),
        ),
        migrations.DeleteModel(
            name='StructureState',
        ),
    ]
