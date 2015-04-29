# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0001_initial'),
        ('residue', '0001_initial'),
        ('mutation', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='mutation',
            name='protein',
            field=models.ForeignKey(to='protein.Protein'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='raw',
            field=models.ForeignKey(to='mutation.MutationRaw'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='refs',
            field=models.ForeignKey(to='mutation.MutationRefs'),
        ),
        migrations.AddField(
            model_name='mutation',
            name='residue',
            field=models.ForeignKey(to='residue.Residue'),
        ),
    ]
