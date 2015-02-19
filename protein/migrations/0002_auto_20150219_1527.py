# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0001_initial'),
        ('protein', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='proteinlinks',
            name='resource',
        ),
        migrations.DeleteModel(
            name='ProteinLinks',
        ),
        migrations.AddField(
            model_name='protein',
            name='web_link',
            field=models.ManyToManyField(to='common.WebLink'),
            preserve_default=True,
        ),
    ]
