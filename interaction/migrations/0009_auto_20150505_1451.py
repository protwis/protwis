# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0008_auto_20150505_1056'),
    ]

    operations = [
        migrations.AlterField(
            model_name='structureligandinteraction',
            name='pdb_file',
            field=models.ForeignKey(to='structure.PdbData', null=True),
        ),
    ]
