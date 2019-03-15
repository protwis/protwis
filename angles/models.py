from django.db import models

class ResidueAngle(models.Model):
    residue     = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    structure   = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    a_angle     = models.FloatField()
    b_angle     = models.FloatField()
    hse         = models.IntegerField(default=0)
    sasa        = models.FloatField(default=0)
    phi         = models.FloatField(default=0, null=True)
    psi         = models.FloatField(default=0, null=True)
    theta       = models.FloatField(default=0, null=True)
    tau         = models.FloatField(default=0, null=True)

    class Meta():
        db_table = 'residue_angles'
        unique_together = ("residue", "structure")
