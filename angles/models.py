from django.db import models

class Angle(models.Model):
    # linked onto the Xtal ProteinConformation, which is linked to the Xtal protein
    residue   = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    angle    = models.FloatField()
    hse      = models.IntegerField(default=0)
    sasa     = models.FloatField(default=0)

    def __str__(self):
        return "angle"

    class Meta:
        unique_together = ("residue", "structure")
