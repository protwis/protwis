# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('protein', '0001_initial'),
        ('structure', '0001_initial'),
        ('ligand', '0001_initial'),
        ('interaction', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='structureligandinteraction',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='interaction_type',
            field=models.ForeignKey(to='interaction.ResidueFragmentInteractionType'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='proteinligandinteraction',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand'),
        ),
        migrations.AddField(
            model_name='proteinligandinteraction',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
        ),
    ]
