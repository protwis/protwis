# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0004_auto_20150911_1726'),
    ]

    operations = [
        migrations.AlterField(
            model_name='proteinfamily',
            name='slug',
            field=models.SlugField(unique=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='proteinsegment',
            name='slug',
            field=models.SlugField(unique=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='proteinsequencetype',
            name='slug',
            field=models.SlugField(unique=True, max_length=20),
        ),
        migrations.AlterField(
            model_name='proteinset',
            name='name',
            field=models.CharField(unique=True, max_length=50),
        ),
        migrations.AlterField(
            model_name='proteinsource',
            name='name',
            field=models.CharField(unique=True, max_length=20),
        ),
        migrations.AlterField(
            model_name='proteinstate',
            name='slug',
            field=models.SlugField(unique=True, max_length=20),
        ),
        migrations.AlterField(
            model_name='species',
            name='latin_name',
            field=models.CharField(unique=True, max_length=100),
        ),
    ]
