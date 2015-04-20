# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0007_auto_20150417_0948'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='proteinfamily',
            options={'ordering': ['id']},
        ),
    ]
