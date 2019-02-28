from django import template

import re

register = template.Library()

@register.filter
def zscales_color_scale ( objs ):
    if objs[2] <= 1:
        return '0'
    else:
        scaling = 10 # manual, should ideally be based on number of observations
        if objs[1]*scaling > 10:
            return '0'
        else:
            return str(int(10 - objs[1]*scaling))
