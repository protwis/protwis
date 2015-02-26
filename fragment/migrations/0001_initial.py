# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0001_initial'),
        ('residue', '0002_residuegenericnumber_protein_segment'),
        ('ligand', '0002_auto_20150225_1441'),
        ('protein', '0003_delete_proteinresource'),
    ]

    operations = [
        migrations.CreateModel(
            name='Fragment',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
            ],
            options={
                'db_table': 'fragment',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Interaction',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('slug', models.CharField(max_length=10)),
                ('description', models.CharField(max_length=200)),
            ],
            options={
                'db_table': 'interaction',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='fragment',
            name='interaction',
            field=models.ForeignKey(to='fragment.Interaction', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='fragment',
            name='ligand',
            field=models.ForeignKey(to='ligand.Ligand', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='fragment',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='fragment',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='fragment',
            name='structure',
            field=models.ForeignKey(to='structure.Structure', null=True),
            preserve_default=True,
        ),
    ]
