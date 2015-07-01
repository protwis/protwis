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
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
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
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('slug', models.CharField(max_length=30, null=True)),
                ('name', models.TextField()),
            ],
            options={
                'db_table': 'publication_journal',
            },
        ),
        migrations.CreateModel(
            name='WebLink',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('index', models.TextField()),
            ],
            options={
                'db_table': 'web_link',
            },
        ),
        migrations.CreateModel(
            name='WebResource',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
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
            field=models.ForeignKey(to='common.WebLink', null=True),
        ),
    ]
