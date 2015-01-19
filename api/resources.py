from tastypie.resources import ModelResource

from protein.models import Protein
from residue.models import Residue


class GetProteinResource(ModelResource):
    class Meta:
        queryset = Protein.objects.all()
        detail_uri_name = 'entry_name'
        allowed_methods = ['get']


class GetSequenceResource(ModelResource):
    class Meta:
        queryset = Protein.objects.all()
        detail_uri_name = 'entry_name'
        fields = ['sequence']
        allowed_methods = ['get']


class GetResidueResource(ModelResource):
    class Meta:
        queryset = Residue.objects.all()
        allowed_methods = ['get']