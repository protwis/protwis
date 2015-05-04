# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0002_auto_20150428_1149'),
    ]

    operations = [
        migrations.CreateModel(
            name='ResidueFragmentAtom',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, serialize=False, primary_key=True)),
                ('atomtype', models.CharField(max_length=20)),
                ('atomnr', models.SmallIntegerField()),
                ('atomclass', models.CharField(max_length=20)),
                ('residuename', models.CharField(max_length=20)),
                ('chain', models.CharField(max_length=20)),
                ('residuenr', models.SmallIntegerField()),
                ('x', models.DecimalField(decimal_places=3, max_digits=6)),
                ('y', models.DecimalField(decimal_places=3, max_digits=6)),
                ('z', models.DecimalField(decimal_places=3, max_digits=6)),
                ('occupancy', models.DecimalField(decimal_places=2, max_digits=6)),
                ('temperature', models.DecimalField(decimal_places=2, max_digits=6)),
                ('element_name', models.CharField(max_length=20)),
                ('interaction', models.ForeignKey(to='interaction.ResidueFragmentInteraction', null=True)),
                ('structureligandpair', models.ForeignKey(to='interaction.StructureLigandInteraction')),
            ],
            options={
                'db_table': 'interaction_residue_fragment_atoms',
            },
        ),
    ]
