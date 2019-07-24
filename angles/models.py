from django.db import models

class ResidueAngle(models.Model):
    residue       = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    structure     = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    a_angle       = models.FloatField(default=0, null=True)
    b_angle       = models.FloatField(default=0, null=True)
    outer_angle   = models.FloatField(default=0, null=True)
    hse           = models.IntegerField(default=0, null=True)
    sasa          = models.FloatField(default=0, null=True)
    rsa           = models.FloatField(default=0, null=True)
    phi           = models.FloatField(default=0, null=True)
    psi           = models.FloatField(default=0, null=True)
    theta         = models.FloatField(default=0, null=True)
    tau           = models.FloatField(default=0, null=True)
    core_distance = models.FloatField(default=0, null=True)
    ss_dssp       = models.CharField(max_length=1, null=True)
    ss_stride     = models.CharField(max_length=1, null=True)

    class Meta():
        db_table = 'residue_angles'
        unique_together = ("residue", "structure")
