# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0003_auto_20150416_1612'),
    ]

    operations = [
        migrations.CreateModel(
            name='ResidueFragmentAtom',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('atomnr', models.SmallIntegerField()),
                ('atomtype', models.CharField(max_length=20)),
                ('residuename', models.CharField(max_length=20)),
                ('chain', models.CharField(max_length=20)),
                ('residuenr', models.SmallIntegerField()),
                ('x', models.DecimalField(decimal_places=3, max_digits=6)),
                ('y', models.DecimalField(decimal_places=3, max_digits=6)),
                ('z', models.DecimalField(decimal_places=3, max_digits=6)),
                ('occupancy', models.DecimalField(decimal_places=3, max_digits=6)),
                ('temperature', models.DecimalField(decimal_places=3, max_digits=6)),
                ('element_name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'interaction_residue_fragment_atoms',
            },
        ),
    ]
