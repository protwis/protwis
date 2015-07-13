from django.shortcuts import render
from rest_framework import views, generics, viewsets
from rest_framework.response import Response
from django.template.loader import render_to_string
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinFamily, Species, ProteinSegment
from residue.models import Residue
from structure.models import Structure
from api.serializers import (ProteinSerializer, ProteinFamilySerializer, SpeciesSerializer, ResidueSerializer,
    ResidueExtendedSerializer, StructureSerializer)
from common.alignment import Alignment

import json

# FIXME add
# similiarity search
# getMutations
# numberPDBfile

def index(request):
    return render(request, 'api/index.html', {'site_title': settings.SITE_TITLE})

class ProteinDetail(generics.RetrieveAPIView):
    """
    Get a single protein instance by entry name
    \n/protein/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    queryset = Protein.objects.filter(sequence_type__slug="wt")
    serializer_class = ProteinSerializer
    lookup_field = 'entry_name'


class ProteinByAccessionDetail(ProteinDetail):
    """
    Get a single protein instance by accession
    \n/protein/accession/{accession}/
    \n{accession} is a protein identifier from Uniprot, e.g. P07550
    """

    lookup_field = 'accession'


class ProteinFamilyList(generics.ListAPIView):
    """
    Get a list of protein families
    \n/proteinfamily/
    """

    queryset = ProteinFamily.objects.all()
    serializer_class = ProteinFamilySerializer


class ProteinFamilyDetail(generics.RetrieveAPIView):
    """
    Get a single protein family instance
    \n/proteinfamily/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    queryset = ProteinFamily.objects.all()
    serializer_class = ProteinFamilySerializer
    lookup_field = 'slug'


class ProteinFamilyChildrenList(generics.ListAPIView):
    """
    Get a list of child families of a protein family
    \n/proteinfamily/children/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    serializer_class = ProteinFamilySerializer
    
    def get_queryset(self):
        family = self.kwargs.get('slug')
        queryset = ProteinFamily.objects.all()
        return queryset.filter(Q(slug__startswith=family) & ~Q(slug=family))


class ProteinsInFamilyList(generics.ListAPIView):
    """
    Get a list of proteins in a protein family
    \n/proteinfamily/proteins/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    serializer_class = ProteinSerializer
    
    def get_queryset(self):
        queryset = Protein.objects.all()
        return queryset.filter(sequence_type__slug='wt', family__slug__startswith=self.kwargs.get('slug'))


class ResiduesList(generics.ListAPIView):
    """
    Get a list of residues of a protein
    \n/residues/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    serializer_class = ResidueSerializer
    
    def get_queryset(self):
        queryset = Residue.objects.all()
        return queryset.filter(protein_conformation__protein__sequence_type__slug='wt',
            protein_conformation__protein__entry_name=self.kwargs.get('entry_name'))


class ResiduesExtendedList(ResiduesList):
    """
    Get a list of residues of a protein, including alternative generic numbers
    \n/residues/extended/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    serializer_class = ResidueExtendedSerializer


class SpeciesList(generics.ListAPIView):
    """
    Get a list of species
    \n/species/
    """

    queryset = Species.objects.all()
    serializer_class = SpeciesSerializer


class SpeciesDetail(generics.RetrieveAPIView):
    """
    Get a single species instance
    \n/species/{latin_name}/
    \n{latin_name} is a species identifier from Uniprot, e.g. Homo Sapiens
    """

    queryset = Species.objects.all()
    serializer_class = SpeciesSerializer
    lookup_field = 'latin_name'


class NumberPDBStructureView(views.APIView):
    """
    WRITEME
    """

    pass


class StructureList(views.APIView):
    """
    Get a list of structures
    \n/structure/
    """

    def get(self, request, pdb_code=None, representative=None):
        structures = self.get_structures(pdb_code, representative)

        # convert objects to a list of dictionaries
        # normal serializers can not be used because of abstraction of tables (e.g. protein_conformation)
        s = []
        for structure in structures:
            s.append({
                'pdb_code': structure.pdb_code.index,
                'protein': structure.protein_conformation.protein.parent.entry_name,
                'preferred_chain': structure.preferred_chain,
                'resolution': structure.resolution,
                'publication_date': structure.publication_date,
                'publication': structure.publication.web_link.__str__(),
                'state': structure.state.name,
                'type': structure.structure_type.name,
            })

        # if a structure is selected, return a single dict rather then a list of dicts
        if len(s) == 1:
            s = s[0]

        return Response(s)

    def get_structures(self, pdb_code=None, representative=None):
        return Structure.objects.all()


class RepresentativeStructureList(StructureList):
    """
    Get a list of representative structures (one for each protein and activation state)
    \n/structure/representative/
    """
    def get_structures(self, pdb_code=None, representative=None):
        return Structure.objects.order_by('protein_conformation__protein__parent', 'state',
            'resolution').distinct('protein_conformation__protein__parent', 'state')

class StructureDetail(StructureList):
    """
    Get a single structure instance
    \n/structure/{pdb_code}/
    \n{pdb_code} is a structure identifier from the Protein Data Bank, e.g. 2RH1
    """
    def get_structures(self, pdb_code=None, representative=None):
        return Structure.objects.filter(pdb_code__index=pdb_code)

class FamilyAlignment(views.APIView):
    """
    Get a full sequence alignment of a protein family
    \n/alignment/family/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    def get(self, request, family=None, segments=None):
        if family is not None:
            ps = Protein.objects.filter(sequence_type__slug='wt', family__slug__startswith=family)
            
            if segments is not None:
                segment_list = segments.split(",")
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()
            
            # render the fasta template as string
            response = render_to_string('alignment/alignment_fasta.html', {'a': a}).split("\n")

            # convert the list to a dict
            ali_dict = {}
            for row in response:
                if row.startswith(">"):
                    k = row[1:]
                else:
                    ali_dict[k] = row
                    k = False

            return Response(ali_dict)


class FamilyAlignmentPartial(FamilyAlignment):
    """
    Get a partial sequence alignment of a protein family
    \n/alignment/family/{slug}/{segments}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{segments} is a comma separated list of protein segment identifiers, e.g. TM3,TM5,TM6
    """


class ProteinAlignment(views.APIView):
    """
    Get a full sequence alignment of two or more proteins
    \n/alignment/protein/{proteins}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human
    """

    def get(self, request, proteins=None, segments=None):
        if proteins is not None:
            protein_list = proteins.split(",")
            ps = Protein.objects.filter(sequence_type__slug='wt', entry_name__in=protein_list)
            
            if segments is not None:
                segment_list = segments.split(",")
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()
            
            # render the fasta template as string
            response = render_to_string('alignment/alignment_fasta.html', {'a': a}).split("\n")

            # convert the list to a dict
            ali_dict = {}
            k = False
            for row in response:
                if row.startswith(">"):
                    k = row[1:]
                elif k:
                    ali_dict[k] = row
                    k = False

            return Response(ali_dict)

class ProteinAlignmentPartial(ProteinAlignment):
    """
    Get a partial sequence alignment of two or more proteins
    \n/alignment/protein/{proteins}/{segments}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human
    \n{segments} is a comma separated list of protein segment identifiers, e.g. TM3,TM5,TM6
    """


class StructureTemplate(views.APIView):
    """
    Get the most similar structure template for a protein using a 7TM alignment
    \n/structure/template/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    def get(self, request, entry_name=None, segments=None):
        if entry_name is not None:
            ref = Protein.objects.get(sequence_type__slug='wt', entry_name=entry_name)

            structures =  Structure.objects.order_by('protein_conformation__protein__parent', 'state',
                'resolution').distinct('protein_conformation__protein__parent', 'state')

            ps = []
            for structure in structures:
                ps.append(structure.protein_conformation.protein.parent)
            
            if segments is not None:
                segment_list = segments.split(",")
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False, category='helix')

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_reference_protein(ref)
            a.load_proteins(ps)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            # calculate identity and similarity of each row compared to the reference
            a.calculate_similarity()
            
            # return the entry_name of the closest template
            return Response(a.proteins[1].protein.entry_name)


class StructureTemplatePartial(StructureTemplate):
    """
    Get the most similar structure template for a protein using a partial alignment
    \n/structure/template/{entry_name}/{segments}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    \n{segments} is a comma separated list of protein segment identifiers, e.g. TM3,TM5,TM6
    """

