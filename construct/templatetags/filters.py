from django import template

register = template.Library()

@register.filter(name='splitGN',is_safe=True)
def splitGN(gn):
    split = gn.split("x")
    with_break = "{}x<br>{}".format(split[0],split[1])
    return with_break

    
@register.filter
def cut_at_20 ( objs ):
    return objs[:20]