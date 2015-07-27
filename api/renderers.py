from rest_framework import renderers

class PDBRenderer(renderers.BaseRenderer):
    media_type = 'chemical/x-pdb'
    format = 'pdb'
    render_style = 'binary'
    filename = 'output.pdb'

    def render(self, data, media_type=None, renderer_context=None):
        return data