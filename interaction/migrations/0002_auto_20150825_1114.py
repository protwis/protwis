# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0001_initial'),
        ('structure', '0001_initial'),
        ('ligand', '0001_initial'),
        ('protein', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='structureligandinteraction',
            name='pdb_file',
            field=models.ForeignKey(null=True, to='structure.PdbData'),
        ),
        migrations.AddField(
            model_name='structureligandinteraction',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='fragment',
            field=models.ForeignKey(to='structure.Fragment'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='interaction_type',
            field=models.ForeignKey(to='interaction.ResidueFragmentInteractionType'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='rotamer',
            field=models.ForeignKey(to='structure.Rotamer'),
        ),
        migrations.AddField(
            model_name='residuefragmentinteraction',
            name='structure_ligand_pair',
            field=models.ForeignKey(to='interaction.StructureLigandInteraction'),
        ),
        migrations.AddField(
            model_name='residuefragmentatom',
            name='interaction',
            field=models.ForeignKey(null=True, to='interaction.ResidueFragmentInteraction'),
        ),
        migrations.AddField(
            model_name='residuefragmentatom',
            name='structureligandpair',
            field=models.ForeignKey(to='interaction.StructureLigandInteraction'),
        ),
        migrations.AddField(
            model_name='proteinligandinteraction',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand'),
        ),
        migrations.AddField(
            model_name='proteinligandinteraction',
            name='protein',
            field=models.ForeignKey(to='protein.ProteinConformation'),
        ),
    ]
