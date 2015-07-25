from django.shortcuts import render

from common.views import AbsTargetSelection


class TargetSelection(AbsTargetSelection):
    pass

class ResidueTablesSelection(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 2
    
    description = 'Select receptors to index by searching or browsing in the middle column. You can select entire receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column, where you can edit the list.\n\nSelect which numbering schemes to use in the middle column. By default, only the GPCRDB numbering scheme is selected.\n\nOnce you have selected all your receptors, click the green button.'


    # Middle section
    numbering_schemes = True

