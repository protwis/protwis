# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0002_structuresegment'),
        ('protein', '0002_auto_20150428_1149'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProteinConformationTemplateStructure',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
            ],
            options={
                'db_table': 'protein_conformation_template_structure',
            },
        ),
        migrations.AddField(
            model_name='proteinconformation',
            name='template_structure',
            field=models.ForeignKey(null=True, to='structure.Structure'),
        ),
        migrations.AddField(
            model_name='proteinconformationtemplatestructure',
            name='protein_conformation',
            field=models.ForeignKey(to='protein.ProteinConformation'),
        ),
        migrations.AddField(
            model_name='proteinconformationtemplatestructure',
            name='protein_segment',
            field=models.ForeignKey(to='protein.ProteinSegment'),
        ),
        migrations.AddField(
            model_name='proteinconformationtemplatestructure',
            name='structure',
            field=models.ForeignKey(to='structure.Structure'),
        ),
    ]
