from django import template

import re

register = template.Library()

@register.filter
def join_attr(obj_list, attr_name, sep=', '):
    return sep.join(getattr(i, attr_name) for i in obj_list)

@register.filter
def only_gproteins ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*G.*", element) and not re.match(".*thase.*|PGS", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def only_arrestins ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*rrestin.*", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def only_fusions ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if not re.match(".*bod.*|.*Ab.*|.*Sign.*|.*G.*|.*restin.*", element) or re.match(".*thase.*|PGS", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def only_antibodies ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*bod.*|.*Ab.*", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def senior_author ( objs ):
    if objs:
        return objs.split(',')[-1]
    else:
        return ''

@register.filter
def cut_at_20 ( objs ):
    return objs[:20]

@register.filter
def get_refined_model_version ( objs ):
    return objs.pdb_data.pdb.split('\n')[0][-10:]

@register.filter
def cut_refined ( objs ):
    return objs.split('_')[0]

@register.filter
def cut_classname ( objs ):
    return objs[5:]