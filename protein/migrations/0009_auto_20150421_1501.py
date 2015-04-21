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
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('name', models.CharField(max_length=100, unique=True)),
                ('sequence', models.TextField(null=True)),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
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
            model_name='proteinfusion',
            name='segment_after',
            field=models.ForeignKey(to='protein.ProteinSegment', related_name='segment_after'),
        ),
        migrations.AddField(
            model_name='proteinfusion',
            name='segment_before',
            field=models.ForeignKey(to='protein.ProteinSegment', related_name='segment_before'),
        ),
    ]
