# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0004_auto_20150327_1120'),
    ]

    operations = [
        migrations.AlterField(
            model_name='residue',
            name='alternative_generic_number',
            field=models.ManyToManyField(related_name='alternative', to='residue.ResidueGenericNumber'),
        ),
    ]
