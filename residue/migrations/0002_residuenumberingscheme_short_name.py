# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuenumberingscheme',
            name='short_name',
            field=models.CharField(default='gpcrdb', max_length=20),
            preserve_default=False,
        ),
    ]
