# Generated by Django 3.1.6 on 2021-05-04 12:09

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0012_auto_20210308_1819'),
    ]

    operations = [
        migrations.RenameField(
            model_name='experimentassay',
            old_name='assay_measure',
            new_name='measured_effector',
        ),
        migrations.AddField(
            model_name='experimentassay',
            name='measured_biological_process',
            field=models.CharField(max_length=60, null=True),
        ),
    ]
