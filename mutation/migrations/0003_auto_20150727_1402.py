# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0002_auto_20150428_1149'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mutation',
            name='ligand',
            field=models.ForeignKey(to='mutation.MutationLigand', null=True),
        ),
    ]
