import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0006_auto_20241031_1314'),
        ('protein', '0002_auto_20180117_1457'),
    ]

    operations = [
        migrations.AddField(
            model_name='structureligandinteraction',
            name='site',
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                to='protein.Site',
            ),
        ),
        migrations.AddField(
            model_name='structureligandinteraction',
            name='chain_res',
            field=models.CharField(max_length=100, null=True),
        ),
    ]
