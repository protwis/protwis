import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0002_auto_20180117_1457'),
        ('signprot', '0012_delete_signprotinteractions'),
    ]

    operations = [
        migrations.AddField(
            model_name='signprotcomplex',
            name='alpha_backbone',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                related_name='alpha_backbone_complex',
                to='protein.Protein',
            ),
        ),
    ]
