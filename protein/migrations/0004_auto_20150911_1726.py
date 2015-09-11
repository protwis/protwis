# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0003_auto_20150902_1535'),
    ]

    operations = [
        migrations.AlterField(
            model_name='gene',
            name='proteins',
            field=models.ManyToManyField(related_name='genes', to='protein.Protein'),
        ),
    ]
