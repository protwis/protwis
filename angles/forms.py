

from django import forms
import datetime
from django.forms import ModelForm, Form
from django.utils.safestring import mark_safe

from functools import partial
DateInput = partial(forms.DateInput, {'class': 'datepicker'})
#controlled vocabularies

