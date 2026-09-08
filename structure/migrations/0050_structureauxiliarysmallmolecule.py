from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0049_structuremodel_model_type'),
    ]

    operations = [
        migrations.CreateModel(
            name='StructureAuxiliarySmallMolecule',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=10)),
                ('title', models.CharField(blank=True, max_length=200, null=True)),
                ('type', models.CharField(max_length=20)),
                ('function', models.CharField(blank=True, max_length=50, null=True)),
                ('chain', models.CharField(max_length=5)),
                ('residue_seq_id', models.IntegerField()),
                ('structure', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='auxiliary_small_molecules', to='structure.structure')),
            ],
            options={
                'db_table': 'structure_auxiliary_small_molecule',
            },
        ),
    ]
