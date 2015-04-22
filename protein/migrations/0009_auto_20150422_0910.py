# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0008_auto_20150417_1312'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProteinFusion',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('name', models.CharField(unique=True, max_length=100)),
                ('sequence', models.TextField(null=True)),
            ],
            options={
                'db_table': 'protein_fusion',
            },
        ),
        migrations.CreateModel(
            name='ProteinFusionProtein',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('protein', models.ForeignKey(to='protein.Protein')),
                ('protein_fusion', models.ForeignKey(to='protein.ProteinFusion')),
            ],
            options={
                'db_table': 'protein_fusion_protein',
            },
        ),
        migrations.AlterModelOptions(
            name='proteinfamily',
            options={'ordering': ('id',)},
        ),
        migrations.AlterModelOptions(
            name='proteinsegment',
            options={'ordering': ('id',)},
        ),
        migrations.RemoveField(
            model_name='proteinsegment',
            name='position',
        ),
        migrations.AddField(
            model_name='proteinfusionprotein',
            name='segment_after',
            field=models.ForeignKey(to='protein.ProteinSegment', related_name='segment_after'),
        ),
        migrations.AddField(
            model_name='proteinfusionprotein',
            name='segment_before',
            field=models.ForeignKey(to='protein.ProteinSegment', related_name='segment_before'),
        ),
        migrations.AddField(
            model_name='proteinfusion',
            name='protein',
            field=models.ManyToManyField(through='protein.ProteinFusionProtein', to='protein.Protein'),
        ),
    ]
