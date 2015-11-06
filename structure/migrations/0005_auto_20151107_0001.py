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
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('reference', models.CharField(max_length=1)),
                ('anomaly', models.ForeignKey(to='protein.ProteinAnomaly', null=True)),
            ],
            options={
                'db_table': 'structure_model_anomalies',
            },
        ),
        migrations.CreateModel(
            name='StructureModelLoopTemplates',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
            ],
            options={
                'db_table': 'structure_model_loop_templates',
            },
        ),
        migrations.CreateModel(
            name='StructureModelResidues',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('origin', models.CharField(max_length=15)),
            ],
            options={
                'db_table': 'structure_model_residues',
            },
        ),
        migrations.AddField(
            model_name='structuremodel',
            name='main_template',
            field=models.ForeignKey(to='structure.Structure', null=True),
        ),
        migrations.AddField(
            model_name='structuremodel',
            name='pdb',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='structuremodel',
            name='state',
            field=models.ForeignKey(to='protein.ProteinState', null=True),
        ),
        migrations.AlterField(
            model_name='structuremodel',
            name='protein',
            field=models.ForeignKey(to='protein.Protein', null=True),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='homology_model',
            field=models.ForeignKey(to='structure.StructureModel', null=True),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='residue',
            field=models.ForeignKey(to='residue.Residue', null=True),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='rotamer',
            field=models.ForeignKey(to='structure.Rotamer', null=True),
        ),
        migrations.AddField(
            model_name='structuremodelresidues',
            name='template',
            field=models.ForeignKey(to='structure.Structure', null=True),
        ),
        migrations.AddField(
            model_name='structuremodellooptemplates',
            name='homology_model',
            field=models.ForeignKey(to='structure.StructureModel', null=True),
        ),
        migrations.AddField(
            model_name='structuremodellooptemplates',
            name='segment',
            field=models.ForeignKey(to='protein.ProteinSegment', null=True),
        ),
        migrations.AddField(
            model_name='structuremodellooptemplates',
            name='template',
            field=models.ForeignKey(to='structure.Structure', null=True),
        ),
        migrations.AddField(
            model_name='structuremodelanomalies',
            name='homology_model',
            field=models.ForeignKey(to='structure.StructureModel', null=True),
        ),
        migrations.AddField(
            model_name='structuremodelanomalies',
            name='template',
            field=models.ForeignKey(to='structure.Structure', null=True),
        ),
    ]
