# Generated by Django 3.0.8 on 2021-01-18 08:18

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0024_auto_20210118_0912'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='analyzedexperiment',
            name='auxiliary_protein',
        ),
        migrations.RemoveField(
            model_name='biasedexperiment',
            name='auxiliary_protein',
        ),
        migrations.AddField(
            model_name='analyzedexperiment',
            name='mutation',
            field=models.CharField(max_length=5, null=True),
        ),
        migrations.AddField(
            model_name='analyzedexperiment',
            name='residue',
            field=models.CharField(max_length=5, null=True),
        ),
    ]
