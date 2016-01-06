# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0004_mutationraw_exp_fold_change'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mutationraw',
            name='exp_mu_effect_value',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationraw',
            name='exp_wt_value',
            field=models.FloatField(),
        ),
    ]
