from django import forms

class PDBform(forms.Form):
    pdbname = forms.CharField(max_length=10)