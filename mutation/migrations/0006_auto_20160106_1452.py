# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0005_auto_20160106_1434'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mutationexperiment',
            name='mu_value',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationexperiment',
            name='wt_value',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationoptional',
            name='mu',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationoptional',
            name='percentage',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationoptional',
            name='wt',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationraw',
            name='exp_fold_change',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationraw',
            name='opt_mu',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationraw',
            name='opt_percentage',
            field=models.FloatField(),
        ),
        migrations.AlterField(
            model_name='mutationraw',
            name='opt_wt',
            field=models.FloatField(),
        ),
    ]
