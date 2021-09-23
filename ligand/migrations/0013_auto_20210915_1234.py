# Generated by Django 3.1.6 on 2021-09-15 10:34

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0005_auto_20210725_1110'),
        ('protein', '0015_proteincouplings_physiological_ligand'),
        ('ligand', '0012_auto_20210308_1819'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='analyzedassay',
            name='t_coefficient',
        ),
        migrations.RemoveField(
            model_name='analyzedassay',
            name='t_value',
        ),
        migrations.RemoveField(
            model_name='analyzedexperiment',
            name='primary',
        ),
        migrations.RemoveField(
            model_name='analyzedexperiment',
            name='secondary',
        ),
        migrations.RemoveField(
            model_name='biasedpathways',
            name='ligand_pubchem',
        ),
        migrations.RemoveField(
            model_name='ligandreceptorstatistics',
            name='ligand',
        ),
        migrations.RemoveField(
            model_name='ligandreceptorstatistics',
            name='reference_protein',
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='delta_emax_ec50',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='delta_relative_transduction_coef',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='effector_family',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='measured_biological_process',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='measured_effector',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='molecule_1',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='molecule_2',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='pathway_level',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='relative_transduction_coef',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='signal_detection_tecnique',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='transduction_coef',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AddField(
            model_name='analyzedexperiment',
            name='external_ligand_ids',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='analyzedexperiment',
            name='receptor_gtpo',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='analyzedexperiment',
            name='receptor_isoform',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='biasedexperiment',
            name='receptor_gtpo',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='biasedexperiment',
            name='receptor_isoform',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='ligandreceptorstatistics',
            name='primary',
            field=models.CharField(max_length=100, null=True),
        ),
        migrations.AddField(
            model_name='ligandreceptorstatistics',
            name='secondary',
            field=models.CharField(max_length=100, null=True),
        ),
        migrations.AlterField(
            model_name='analyzedassay',
            name='potency',
            field=models.FloatField(max_length=60, null=True),
        ),
        migrations.AlterField(
            model_name='analyzedassay',
            name='signalling_protein',
            field=models.CharField(max_length=60, null=True),
        ),
        migrations.AlterField(
            model_name='biasedexperiment',
            name='auxiliary_protein',
            field=models.TextField(null=True),
        ),
        migrations.CreateModel(
            name='GTP_endogenous_ligand',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ligand_type', models.TextField(null=True)),
                ('endogenous_princip', models.TextField(null=True)),
                ('pec50_avg', models.FloatField(max_length=60, null=True)),
                ('pec50_min', models.FloatField(max_length=60, null=True)),
                ('pec50_max', models.FloatField(max_length=60, null=True)),
                ('pKi_avg', models.FloatField(max_length=60, null=True)),
                ('pKi_min', models.FloatField(max_length=60, null=True)),
                ('pKi_max', models.FloatField(max_length=60, null=True)),
                ('gpt_link', models.TextField(null=True)),
                ('ligand', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.ligand')),
                ('publication', models.ManyToManyField(null=True, to='common.Publication')),
                ('receptor', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.protein')),
            ],
        ),
        migrations.CreateModel(
            name='BiasedExperimentAssay',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('signalling_protein', models.CharField(max_length=60, null=True)),
                ('family', models.CharField(max_length=60, null=True)),
                ('cell_line', models.CharField(max_length=60, null=True)),
                ('assay_type', models.CharField(max_length=60, null=True)),
                ('molecule_1', models.CharField(max_length=60, null=True)),
                ('molecule_2', models.CharField(max_length=60, null=True)),
                ('pathway_level', models.TextField(null=True)),
                ('measured_biological_process', models.CharField(max_length=60, null=True)),
                ('signal_detection_tecnique', models.TextField(null=True)),
                ('assay_time_resolved', models.CharField(max_length=60, null=True)),
                ('ligand_function', models.CharField(max_length=60, null=True)),
                ('quantitive_measure_type', models.CharField(max_length=60, null=True)),
                ('quantitive_activity', models.FloatField(max_length=60, null=True)),
                ('quantitive_activity_sign', models.CharField(max_length=3, null=True)),
                ('quantitive_unit', models.CharField(max_length=60, null=True)),
                ('qualitative_activity', models.CharField(max_length=60, null=True)),
                ('quantitive_efficacy', models.FloatField(max_length=60, null=True)),
                ('efficacy_measure_type', models.CharField(max_length=60, null=True)),
                ('efficacy_sign', models.CharField(max_length=3, null=True)),
                ('efficacy_unit', models.CharField(max_length=60, null=True)),
                ('bias_reference', models.CharField(max_length=60, null=True)),
                ('transduction_coef', models.FloatField(max_length=60, null=True)),
                ('relative_transduction_coef', models.FloatField(max_length=60, null=True)),
                ('biased_experiment', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='experiment_data', to='ligand.biasedexperiment')),
                ('emax_ligand_reference', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='ExperimentAssay.bias_ligand_reference+', to='ligand.ligand')),
            ],
        ),
        migrations.AddField(
            model_name='analyzedassay',
            name='reference_assay_initial',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='test_ExperimentAssay.bias_ligand_reference_assay+', to='ligand.biasedexperimentassay'),
        ),
        migrations.AlterField(
            model_name='experimentassayauthors',
            name='experiment',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='experiment_data_authors', to='ligand.biasedexperimentassay'),
        ),
        migrations.DeleteModel(
            name='ExperimentAssay',
        ),
    ]
