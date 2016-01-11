# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0008_proteinsegment_fully_aligned'),
        ('structure', '0005_auto_20151113_1551'),
    ]

    operations = [
        migrations.CreateModel(
            name='StructureCoordinates',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'structure_coordinates',
            },
        ),
        migrations.CreateModel(
            name='StructureCoordinatesDescription',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('text', models.TextField()),
            ],
            options={
                'db_table': 'structure_coordinates_description',
            },
        ),
        migrations.CreateModel(
            name='StructureEngineering',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'structure_engineering',
            },
        ),
        migrations.CreateModel(
            name='StructureEngineeringDescription',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('text', models.TextField()),
            ],
            options={
                'db_table': 'structure_engineering_description',
            },
        ),
        migrations.CreateModel(
            name='StructureSegmentModeling',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('start', models.IntegerField()),
                ('end', models.IntegerField()),
                ('protein_segment', models.ForeignKey(to='protein.ProteinSegment')),
                ('structure', models.ForeignKey(to='structure.Structure')),
            ],
            options={
                'db_table': 'structure_segment_modeling',
            },
        ),
        migrations.AddField(
            model_name='structureengineering',
            name='description',
            field=models.ForeignKey(to='structure.StructureEngineeringDescription'),
        ),
        migrations.AddField(
            model_name='structureengineering',
            name='protein_segment',
            field=models.ForeignKey(to='protein.ProteinSegment'),
        ),
        migrations.AddField(
            model_name='structureengineering',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='structurecoordinates',
            name='description',
            field=models.ForeignKey(to='structure.StructureCoordinatesDescription'),
        ),
        migrations.AddField(
            model_name='structurecoordinates',
            name='protein_segment',
            field=models.ForeignKey(to='protein.ProteinSegment'),
        ),
        migrations.AddField(
            model_name='structurecoordinates',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
    ]
