# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0006_auto_20160111_1008'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structurecoordinatesdescription',
            name='text',
            field=models.CharField(unique=True, max_length=200),
        ),
        migrations.AlterField(
            model_name='structureengineeringdescription',
            name='text',
            field=models.CharField(unique=True, max_length=200),
        ),
    ]
