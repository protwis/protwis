from django.db import models

class Angle(models.Model):
    residue   = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    angle    = models.FloatField()
    b_angle  = models.FloatField()
    diff_med = models.FloatField()
    sign_med = models.FloatField(default=0)
    hse      = models.IntegerField(default=0)
    sasa     = models.FloatField(default=0)

    def __str__(self):
        return "angle"

    class Meta:
        unique_together = ("residue", "structure")
