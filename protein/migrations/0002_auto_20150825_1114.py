# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('structure', '0001_initial'),
        ('ligand', '0001_initial'),
        ('common', '0001_initial'),
        ('protein', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='proteinconformationtemplatestructure',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='protein_anomalies',
            field=models.ManyToManyField(to='protein.ProteinAnomaly'),
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='state',
            field=models.ForeignKey(to='protein.ProteinState'),
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='template_structure',
            field=models.ForeignKey(null=True, to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='proteinanomalyruleset',
            name='protein_anomaly',
            field=models.ForeignKey(related_name='rulesets', to='protein.ProteinAnomaly'),
        ),
        migrations.AddField(
            model_name='proteinanomalyrule',
            name='generic_number',
            field=models.ForeignKey(to='residue.ResidueGenericNumber'),
        ),
        migrations.AddField(
            model_name='proteinanomalyrule',
            name='rule_set',
            field=models.ForeignKey(related_name='rules', to='protein.ProteinAnomalyRuleSet'),
        ),
        migrations.AddField(
            model_name='proteinanomaly',
            name='anomaly_type',
            field=models.ForeignKey(to='protein.ProteinAnomalyType'),
        ),
        migrations.AddField(
            model_name='proteinanomaly',
            name='generic_number',
            field=models.ForeignKey(to='residue.ResidueGenericNumber'),
        ),
        migrations.AddField(
            model_name='proteinalias',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='protein',
            name='endogenous_ligands',
            field=models.ManyToManyField(to='ligand.Ligand'),
        ),
        migrations.AddField(
            model_name='protein',
            name='family',
            field=models.ForeignKey(to='protein.ProteinFamily'),
        ),
        migrations.AddField(
            model_name='protein',
            name='parent',
            field=models.ForeignKey(null=True, to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='protein',
            name='residue_numbering_scheme',
            field=models.ForeignKey(to='residue.ResidueNumberingScheme'),
        ),
        migrations.AddField(
            model_name='protein',
            name='sequence_type',
            field=models.ForeignKey(to='protein.ProteinSequenceType'),
        ),
        migrations.AddField(
            model_name='protein',
            name='source',
            field=models.ForeignKey(to='protein.ProteinSource'),
        ),
        migrations.AddField(
            model_name='protein',
            name='species',
            field=models.ForeignKey(to='protein.Species'),
        ),
        migrations.AddField(
            model_name='protein',
            name='states',
            field=models.ManyToManyField(to='protein.ProteinState', through='protein.ProteinConformation'),
        ),
        migrations.AddField(
            model_name='protein',
            name='web_links',
            field=models.ManyToManyField(to='common.WebLink'),
        ),
        migrations.AddField(
            model_name='gene',
            name='proteins',
            field=models.ManyToManyField(to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='gene',
            name='species',
            field=models.ForeignKey(to='protein.Species'),
        ),
    ]
