# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Publication',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('title', models.TextField()),
                ('year', models.IntegerField()),
                ('citation', models.TextField()),
            ],
            options={
                'db_table': 'publication',
            },
        ),
        migrations.CreateModel(
            name='PublicationJournal',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('slug', models.CharField(max_length=30, null=True)),
                ('name', models.TextField()),
            ],
            options={
                'db_table': 'publication_journal',
            },
        ),
        migrations.AddField(
            model_name='publication',
            name='journal',
            field=models.ForeignKey(to='common.PublicationJournal'),
        ),
        migrations.AddField(
            model_name='publication',
            name='web_link',
            field=models.ForeignKey(to='common.WebLink'),
        ),
    ]
