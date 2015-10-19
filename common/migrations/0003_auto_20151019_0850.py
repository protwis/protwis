# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0002_auto_20150927_1705'),
    ]

    operations = [
        migrations.AlterField(
            model_name='publicationjournal',
            name='slug',
            field=models.CharField(max_length=200, null=True),
        ),
        migrations.AlterUniqueTogether(
            name='weblink',
            unique_together=set([('web_resource', 'index')]),
        ),
    ]
