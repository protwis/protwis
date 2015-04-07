# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('protein', '0001_initial'),
        ('common', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='protein',
            name='residue_numbering_scheme',
            field=models.ForeignKey(to='residue.ResidueNumberingScheme'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='protein',
            name='source',
            field=models.ForeignKey(to='protein.ProteinSource'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='protein',
            name='species',
            field=models.ForeignKey(to='protein.Species'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='protein',
            name='web_link',
            field=models.ManyToManyField(to='common.WebLink'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='gene',
            name='proteins',
            field=models.ManyToManyField(to='protein.Protein'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='gene',
            name='species',
            field=models.ForeignKey(to='protein.Species'),
            preserve_default=True,
        ),
    ]
