# Generated by Django 3.0.3 on 2024-12-04 09:47

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0026_smiles_search_rdkit'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ligand',
            name='pdbe',
            field=models.CharField(max_length=5, null=True),
        ),
    ]