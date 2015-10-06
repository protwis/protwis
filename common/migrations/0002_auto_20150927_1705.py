# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='ReleaseNotes',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('date', models.DateField()),
                ('html', models.TextField()),
            ],
            options={
                'ordering': ('-date',),
                'db_table': 'release_notes',
            },
        ),
        migrations.CreateModel(
            name='ReleaseStatistics',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('value', models.IntegerField()),
                ('release', models.ForeignKey(to='common.ReleaseNotes')),
            ],
            options={
                'ordering': ('id',),
                'db_table': 'release_statistics',
            },
        ),
        migrations.CreateModel(
            name='ReleaseStatisticsType',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('name', models.CharField(max_length=200)),
            ],
            options={
                'db_table': 'release_statistics_type',
            },
        ),
        migrations.AddField(
            model_name='releasestatistics',
            name='statistics_type',
            field=models.ForeignKey(to='common.ReleaseStatisticsType'),
        ),
    ]
