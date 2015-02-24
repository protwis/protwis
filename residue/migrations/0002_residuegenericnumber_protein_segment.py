# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        #('protein', '0003_delete_proteinresource'),
        ('residue', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='residuegenericnumber',
            name='protein_segment',
            field=models.ForeignKey(to='protein.ProteinSegment', null=True),
            preserve_default=True,
        ),
    ]
