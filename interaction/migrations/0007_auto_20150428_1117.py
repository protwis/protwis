# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0006_residuefragmentatom_atomclass'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuefragmentatom',
            name='StructureLigandPair',
            field=models.ForeignKey(default='', to='interaction.StructureLigandInteraction'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='residuefragmentatom',
            name='interaction',
            field=models.ForeignKey(to='interaction.ResidueFragmentInteraction', null=True),
        ),
    ]
