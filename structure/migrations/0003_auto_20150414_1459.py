# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0003_auto_20150412_2031'),
        ('structure', '0002_auto_20150408_2100'),
    ]

    operations = [
        migrations.CreateModel(
            name='StructureModel',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('protein', models.ForeignKey(to='protein.Protein')),
            ],
            options={
                'db_table': 'structure_model',
            },
        ),
        migrations.RenameField(
            model_name='structure',
            old_name='pdb_publication_date',
            new_name='publication_date',
        ),
        migrations.RemoveField(
            model_name='structure',
            name='endogenous_ligand',
        ),
        migrations.RemoveField(
            model_name='structure',
            name='xray_ligand',
        ),
        migrations.RemoveField(
            model_name='structuretype',
            name='description',
        ),
        migrations.AddField(
            model_name='structurestabilizingagent',
            name='name',
            field=models.CharField(default='', max_length=100),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='structuretype',
            name='name',
            field=models.CharField(default='', max_length=100),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='structure',
            name='stabilizing_agents',
            field=models.ManyToManyField(to='structure.StructureStabilizingAgent'),
        ),
    ]
