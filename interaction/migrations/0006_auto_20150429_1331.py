# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0005_remove_residuefragmentinteraction_fragment'),
    ]

    operations = [
        migrations.RenameField(
            model_name='residuefragmentinteraction',
            old_name='sfragment',
            new_name='fragment',
        ),
    ]
