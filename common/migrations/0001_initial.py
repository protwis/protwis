# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Publication',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('title', models.TextField()),
                ('authors', models.TextField()),
                ('year', models.IntegerField()),
                ('reference', models.TextField()),
            ],
            options={
                'db_table': 'publication',
            },
        ),
        migrations.CreateModel(
            name='PublicationJournal',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('slug', models.CharField(null=True, max_length=30)),
                ('name', models.TextField()),
            ],
            options={
                'db_table': 'publication_journal',
            },
        ),
        migrations.CreateModel(
            name='WebLink',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('index', models.TextField()),
            ],
            options={
                'db_table': 'web_link',
            },
        ),
        migrations.CreateModel(
            name='WebResource',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(max_length=200, default='')),
                ('url', models.TextField()),
            ],
            options={
                'db_table': 'web_resource',
            },
        ),
        migrations.AddField(
            model_name='weblink',
            name='web_resource',
            field=models.ForeignKey(to='common.WebResource'),
        ),
        migrations.AddField(
            model_name='publication',
            name='journal',
            field=models.ForeignKey(to='common.PublicationJournal'),
        ),
        migrations.AddField(
            model_name='publication',
            name='web_link',
            field=models.ForeignKey(null=True, to='common.WebLink'),
        ),
    ]
