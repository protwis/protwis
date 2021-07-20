# Generated by Django 3.1.6 on 2021-07-20 13:09

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0014_gtp_endogenous_ligand'),
    ]

    operations = [
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='gpt_link',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='pKi_avg',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='pKi_max',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='pKi_min',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='pec50_avg',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='pec50_max',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='gtp_endogenous_ligand',
            name='pec50_min',
            field=models.FloatField(max_length=60, null=True),
        ),
    ]
