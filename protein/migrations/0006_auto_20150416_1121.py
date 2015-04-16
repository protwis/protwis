# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0005_auto_20150415_2149'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='proteinanomalyruleset',
            name='protein_anomaly_rules',
        ),
        migrations.AddField(
            model_name='proteinanomalyrule',
            name='rule_set',
            field=models.ForeignKey(to='protein.ProteinAnomalyRuleSet', default=0),
            preserve_default=False,
        ),
        migrations.RemoveField(
            model_name='proteinanomalyrule',
            name='amino_acid',
        ),
        migrations.AddField(
            model_name='proteinanomalyrule',
            name='amino_acid',
            field=models.CharField(default='', max_length=1),
            preserve_default=False,
        ),
    ]
