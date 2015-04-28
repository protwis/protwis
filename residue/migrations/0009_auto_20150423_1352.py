# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0010_auto_20150423_1352'),
        ('residue', '0008_auto_20150417_1312'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='residue',
            name='protein',
        ),
        migrations.AddField(
            model_name='residue',
            name='protein_conformation',
            field=models.ForeignKey(to='protein.ProteinConformation', default=0),
            preserve_default=False,
        ),
    ]
