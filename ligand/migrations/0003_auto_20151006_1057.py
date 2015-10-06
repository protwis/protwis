# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0002_auto_20150915_2308'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='ligand',
            unique_together=set([('name', 'canonical')]),
        ),
    ]
