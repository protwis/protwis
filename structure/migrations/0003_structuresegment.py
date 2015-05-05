# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0003_auto_20150505_1302'),
        ('structure', '0002_auto_20150429_1313'),
    ]

    operations = [
        migrations.CreateModel(
            name='StructureSegment',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('start', models.IntegerField()),
                ('end', models.IntegerField()),
                ('protein_segment', models.ForeignKey(to='protein.ProteinSegment')),
                ('structure', models.ForeignKey(to='structure.Structure')),
            ],
            options={
                'db_table': 'structure_segment',
            },
        ),
    ]
