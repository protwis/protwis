# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0007_auto_20150429_1408'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='residuefragmentinteraction',
            name='ligand',
        ),
        migrations.RemoveField(
            model_name='residuefragmentinteraction',
            name='residue',
        ),
        migrations.RemoveField(
            model_name='residuefragmentinteraction',
            name='structure',
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='structure_ligand_pair',
            field=models.ForeignKey(to='interaction.StructureLigandInteraction', default=0),
            preserve_default=False,
        ),
    ]
