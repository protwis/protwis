# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0002_auto_20150429_1313'),
        ('interaction', '0006_auto_20150429_1331'),
    ]

    operations = [
        migrations.AddField(
            model_name='structureligandinteraction',
            name='pdb_file',
            field=models.ForeignKey(to='structure.PdbData', default=0),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='proteinligandinteraction',
            name='protein',
            field=models.ForeignKey(to='protein.ProteinConformation'),
        ),
    ]
