# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0002_auto_20150316_1232'),
    ]

    operations = [
        migrations.CreateModel(
            name='ProteinSequenceType',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('slug', models.SlugField(unique=True, max_length=20)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'protein_sequence_type',
            },
        ),
        migrations.AddField(
            model_name='protein',
            name='parent',
            field=models.ForeignKey(null=True, to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='protein',
            name='sequence_type',
            field=models.ForeignKey(default=1, to='protein.ProteinSequenceType'),
            preserve_default=False,
        ),
    ]
