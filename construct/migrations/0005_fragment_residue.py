# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0003_residuenumberingscheme_parent'),
        ('structure', '0004_merge'),
    ]

    operations = [
        migrations.AddField(
            model_name='fragment',
            name='residue',
            field=models.ForeignKey(default='', to='residue.Residue'),
            preserve_default=False,
        ),
    ]
