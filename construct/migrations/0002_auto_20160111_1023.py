# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('construct', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelTable(
            name='auxprotein',
            table='construct_aux_protein',
        ),
        migrations.AlterModelTable(
            name='auxproteintype',
            table='construct_aux_protein_type',
        ),
        migrations.AlterModelTable(
            name='chemical',
            table='construct_chemical',
        ),
        migrations.AlterModelTable(
            name='chemicalconc',
            table='construct_chemical_conc',
        ),
        migrations.AlterModelTable(
            name='chemicallist',
            table='construct_chemical_list',
        ),
        migrations.AlterModelTable(
            name='chemicalmodification',
            table='construct_chemical_modification',
        ),
        migrations.AlterModelTable(
            name='chemicaltype',
            table='construct_chemical_type',
        ),
        migrations.AlterModelTable(
            name='constructcrystallizationligandconc',
            table='construct_ligand_conc_of_crystallization',
        ),
        migrations.AlterModelTable(
            name='crystallizationmethodtypes',
            table='construct_crystallization_method_types',
        ),
        migrations.AlterModelTable(
            name='purificationstep',
            table='construct_purification_step',
        ),
        migrations.AlterModelTable(
            name='purificationsteptype',
            table='construct_purification_step_type',
        ),
    ]
