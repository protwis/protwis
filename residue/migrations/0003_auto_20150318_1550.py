# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0002_residuenumberingscheme_short_name'),
    ]

    operations = [
        migrations.AlterField(
            model_name='residue',
            name='alternative_generic_number',
            field=models.ManyToManyField(related_name='alternative', to='residue.ResidueGenericNumber', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='residue',
            name='display_generic_number',
            field=models.ForeignKey(related_name='display', to='residue.ResidueGenericNumber', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='residue',
            name='generic_number',
            field=models.ForeignKey(related_name='compare', to='residue.ResidueGenericNumber', null=True),
            preserve_default=True,
        ),
    ]
