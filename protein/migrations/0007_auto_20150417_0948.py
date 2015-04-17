# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0006_auto_20150416_1121'),
    ]

    operations = [
        migrations.AlterField(
            model_name='protein',
            name='accession',
            field=models.CharField(db_index=True, null=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='protein',
            name='entry_name',
            field=models.SlugField(null=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='protein',
            name='source',
            field=models.ForeignKey(to='protein.ProteinSource', null=True),
        ),
        migrations.AlterField(
            model_name='proteinanomalytype',
            name='slug',
            field=models.SlugField(max_length=20),
        ),
        migrations.AlterField(
            model_name='proteinfamily',
            name='slug',
            field=models.SlugField(max_length=100),
        ),
        migrations.AlterField(
            model_name='proteinsegment',
            name='slug',
            field=models.SlugField(max_length=100),
        ),
        migrations.AlterField(
            model_name='proteinsequencetype',
            name='slug',
            field=models.SlugField(max_length=20),
        ),
    ]
