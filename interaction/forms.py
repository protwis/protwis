from django import forms

class PDBform(forms.Form):
    pdbname = forms.CharField(max_length=10, required=False)
    file = forms.FileField(label='Select a file', help_text='max. 42 megabytes', required=False)