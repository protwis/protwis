# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Documentation',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('title', models.TextField()),
                ('description', models.TextField()),
                ('image', models.TextField()),
                ('html', models.TextField()),
            ],
            options={
                'db_table': 'documentation',
            },
        ),
    ]
