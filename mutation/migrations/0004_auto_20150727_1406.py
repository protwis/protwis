# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mutation', '0003_auto_20150727_1402'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mutation',
            name='refs',
            field=models.ForeignKey(null=True, to='mutation.MutationRefs'),
        ),
    ]
