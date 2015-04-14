# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0002_auto_20150413_1126'),
    ]

    operations = [
        migrations.RenameField(
            model_name='publication',
            old_name='citation',
            new_name='reference',
        ),
        migrations.AddField(
            model_name='publication',
            name='authors',
            field=models.TextField(default=''),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='publication',
            name='web_link',
            field=models.ForeignKey(to='common.WebLink', null=True),
        ),
    ]
