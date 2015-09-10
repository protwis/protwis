from rest_framework import views
from rest_framework.response import Response

from api.views import *

import json
import logging


class HelixBoxView(views.APIView):
    """
    Get SVG source code for a protein's helix box plot
    \n/plot/helixbox/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    def get(self, request, entry_name=None):
        if entry_name is not None:
            p = Protein.objects.get(entry_name=entry_name)

            return Response(str(p.get_helical_box()).split('\n'))


class SnakePlotView(views.APIView):
    """
    Get SVG source code for a protein's snake plot
    \n/plot/snake/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    def get(self, request, entry_name=None):
        if entry_name is not None:
            p = Protein.objects.get(entry_name=entry_name)

            return Response(str(p.get_snake_plot()).split('\n'))