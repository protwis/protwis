# -*- coding: utf-8 -*-
# Generated by Django 1.11.2 on 2017-09-08 05:58
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('ligand', '0001_initial'),
        ('interaction', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='structureligandinteraction',
            name='ligand',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.Ligand'),
        ),
        migrations.AddField(
            model_name='structureligandinteraction',
            name='ligand_role',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.LigandRole'),
        ),
    ]
