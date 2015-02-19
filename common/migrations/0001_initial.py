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
                ('id', models.AutoField(serialize=False, verbose_name='ID', primary_key=True, auto_created=True)),
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
                ('id', models.AutoField(serialize=False, verbose_name='ID', primary_key=True, auto_created=True)),
                ('name', models.CharField(max_length=200)),
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
