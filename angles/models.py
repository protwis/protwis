from django.db import models

class ResidueAngle(models.Model):
    residue     = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    structure   = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    angle       = models.FloatField()
    b_angle     = models.FloatField()
    hse         = models.IntegerField(default=0)
    sasa        = models.FloatField(default=0)
    phi         = models.FloatField(default=0)
    psi         = models.FloatField(default=0)
    theta       = models.FloatField(default=0)
    tau         = models.FloatField(default=0)

    class Meta():
        db_table = 'residue_angles'
        unique_together = ("residue", "structure")
