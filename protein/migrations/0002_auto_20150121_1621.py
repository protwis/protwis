# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelTable(
            name='gene',
            table='gene',
        ),
        migrations.AlterModelTable(
            name='protein',
            table='protein',
        ),
        migrations.AlterModelTable(
            name='proteinalias',
            table='protein_alias',
        ),
        migrations.AlterModelTable(
            name='proteinfamily',
            table='protein_family',
        ),
        migrations.AlterModelTable(
            name='proteinlinks',
            table='protein_links',
        ),
        migrations.AlterModelTable(
            name='proteinresource',
            table='protein_reosurce',
        ),
        migrations.AlterModelTable(
            name='proteinsegment',
            table='protein_segment',
        ),
        migrations.AlterModelTable(
            name='proteinset',
            table='protein_set',
        ),
        migrations.AlterModelTable(
            name='proteinsource',
            table='protein_source',
        ),
        migrations.AlterModelTable(
            name='species',
            table='species',
        ),
    ]
