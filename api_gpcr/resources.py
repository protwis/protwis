from tastypie.resources import ModelResource

from protein.models import Species


class GetSpeciesResource(ModelResource):
    class Meta:
        queryset = Species.objects.all()
        allowed_methods = ['get']