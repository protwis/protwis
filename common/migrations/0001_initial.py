# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='WebLink',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('index', models.TextField()),
            ],
            options={
                'db_table': 'web_link',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='WebResource',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('slug', models.SlugField(max_length=20)),
                ('name', models.CharField(default='', max_length=200)),
                ('url', models.TextField()),
            ],
            options={
                'db_table': 'web_resource',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='weblink',
            name='web_resource',
            field=models.ForeignKey(to='common.WebResource'),
            preserve_default=True,
        ),
    ]
