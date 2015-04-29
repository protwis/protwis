# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0004_residuefragmentatom'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuefragmentatom',
            name='interaction',
            field=models.ForeignKey(to='interaction.ResidueFragmentInteraction', default=''),
            preserve_default=False,
        ),
    ]
