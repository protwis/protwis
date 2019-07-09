from django import template
register = template.Library()


@register.filter
def replace_slash(string):
    if string =='pERK1/2 activation':
        string=='pERK1_2 activation'
    else:
        pass
    return string
