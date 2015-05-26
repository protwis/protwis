# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0004_proteinsegment_partial'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='proteinanomalyruleset',
            options={'ordering': ('id',)},
        ),
        migrations.AlterField(
            model_name='proteinanomalyrule',
            name='rule_set',
            field=models.ForeignKey(to='protein.ProteinAnomalyRuleSet', related_name='rules'),
        ),
        migrations.AlterField(
            model_name='proteinanomalyruleset',
            name='protein_anomaly',
            field=models.ForeignKey(to='protein.ProteinAnomaly', related_name='rulesets'),
        ),
    ]
