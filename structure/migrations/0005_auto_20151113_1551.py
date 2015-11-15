# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0008_proteinsegment_fully_aligned'),
        ('residue', '0002_auto_20150902_2223'),
        ('structure', '0004_auto_20151006_1057'),
    ]

    operations = [
        migrations.CreateModel(
            name='StructureModelAnomalies',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('reference', models.CharField(max_length=1)),
                ('anomaly', models.ForeignKey(to='protein.ProteinAnomaly')),
            ],
            options={
                'db_table': 'structure_model_anomalies',
            },
        ),
        migrations.CreateModel(
            name='StructureModelLoopTemplates',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
            ],
            options={
                'db_table': 'structure_model_loop_templates',
            },
        ),
        migrations.CreateModel(
            name='StructureModelResidues',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('sequence_number', models.IntegerField()),
                ('origin', models.CharField(max_length=15)),
            ],
            options={
                'db_table': 'structure_model_residues',
            },
        ),
        migrations.AddField(
            model_name='structuremodel',
            name='main_template',
            field=models.ForeignKey(default=0, to='structure.Structure'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='structuremodel',
            name='pdb',
            field=models.TextField(default=''),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='structuremodel',
            name='state',
            field=models.ForeignKey(default=0, to='protein.ProteinState'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='homology_model',
            field=models.ForeignKey(to='structure.StructureModel'),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='rotamer',
            field=models.ForeignKey(to='structure.Rotamer', null=True),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='segment',
            field=models.ForeignKey(to='protein.ProteinSegment'),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='template',
            field=models.ForeignKey(to='structure.Structure', null=True),
        ),
        migrations.AddField(
            model_name='structuremodellooptemplates',
            name='homology_model',
            field=models.ForeignKey(to='structure.StructureModel'),
        ),
        migrations.AddField(
            model_name='structuremodellooptemplates',
            name='segment',
            field=models.ForeignKey(to='protein.ProteinSegment'),
        ),
        migrations.AddField(
            model_name='structuremodellooptemplates',
            name='template',
            field=models.ForeignKey(to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='structuremodelanomalies',
            name='homology_model',
            field=models.ForeignKey(to='structure.StructureModel'),
        ),
        migrations.AddField(
            model_name='structuremodelanomalies',
            name='template',
            field=models.ForeignKey(to='structure.Structure'),
        ),
    ]
