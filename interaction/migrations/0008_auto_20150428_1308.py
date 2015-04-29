# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0007_auto_20150428_1117'),
    ]

    operations = [
        migrations.RenameField(
            model_name='residuefragmentatom',
            old_name='StructureLigandPair',
            new_name='structureligandpair',
        ),
    ]
