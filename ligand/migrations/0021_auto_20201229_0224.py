# Generated by Django 3.0.8 on 2020-12-29 01:24

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0020_auto_20201130_1437'),
    ]

    operations = [
        migrations.AlterField(
            model_name='assayexperiment',
            name='standard_value',
            field=models.DecimalField(decimal_places=1, max_digits=20, null=True),
        ),
        migrations.CreateModel(
            name='AssayExperimentSource',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('database', models.CharField(max_length=20, null=True)),
                ('database_id', models.CharField(max_length=30, null=True)),
                ('assay', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.AssayExperiment')),
            ],
        ),
    ]
