# Generated by Django 3.1.6 on 2021-02-16 11:01

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0010_ligandreceptorstatistics'),
    ]

    operations = [
        migrations.RenameField(
            model_name='biasedpathways',
            old_name='lignad_pubchem',
            new_name='ligand_pubchem',
        ),
    ]
