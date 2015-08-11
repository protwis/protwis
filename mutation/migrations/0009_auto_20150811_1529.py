# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0008_auto_20150811_1507'),
    ]

    operations = [
        migrations.RenameField(
            model_name='mutation',
            old_name='mutation_to',
            new_name='amino_acid',
        ),
    ]
