# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0003_auto_20150825_1636'),
    ]

    operations = [
        migrations.AddField(
            model_name='mutationraw',
            name='exp_fold_change',
            field=models.DecimalField(max_digits=10, default=0, decimal_places=2),
            preserve_default=False,
        ),
    ]
