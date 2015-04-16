# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0005_structure_protein_anomaly'),
    ]

    operations = [
        migrations.CreateModel(
            name='StructureState',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'structure_state',
            },
        ),
        migrations.RenameField(
            model_name='structure',
            old_name='protein_anomaly',
            new_name='protein_anomalies',
        ),
        migrations.AddField(
            model_name='structure',
            name='state',
            field=models.ForeignKey(to='structure.StructureState', default=0),
            preserve_default=False,
        ),
    ]
