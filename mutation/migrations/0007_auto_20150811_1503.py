# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0006_auto_20150810_1554'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='mutationexperiment',
            name='mutation_to',
        ),
        migrations.AddField(
            model_name='mutationexperiment',
            name='mutation',
            field=models.ForeignKey(default='', to='mutation.Mutation'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='mutation',
            name='mutation_type',
            field=models.ForeignKey(to='mutation.MutationType', null=True),
        ),
    ]
