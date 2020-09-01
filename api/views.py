from django.shortcuts import render
from rest_framework import views, generics, viewsets
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser, FormParser, FileUploadParser
from rest_framework.renderers import JSONRenderer
from django.template.loader import render_to_string
from django.db.models import Q
from django.conf import settings

from interaction.models import ResidueFragmentInteraction
from mutation.models import MutationRaw
from protein.models import Protein, ProteinConformation, ProteinFamily, Species, ProteinSegment
from residue.models import Residue, ResidueGenericNumber, ResidueNumberingScheme, ResidueGenericNumberEquivalent
from structure.models import Structure
from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.sequence_parser import SequenceParser
from api.serializers import (ProteinSerializer, ProteinFamilySerializer, SpeciesSerializer, ResidueSerializer,
                             ResidueExtendedSerializer, StructureSerializer,
                             StructureLigandInteractionSerializer,
                             MutationSerializer)
from api.renderers import PDBRenderer
from common.alignment import Alignment
from common.definitions import *
from drugs.models import Drugs

import json, os
from io import StringIO
from Bio.PDB import PDBIO, parse_pdb_header
from collections import OrderedDict

# FIXME add
# getMutations
# numberPDBfile
import coreapi
from urllib.parse import urlparse
from urllib.parse import urljoin
from rest_framework import renderers, response, schemas
from rest_framework.decorators import api_view, renderer_classes
from rest_framework import response, schemas
from rest_framework_swagger.views import get_swagger_view

schema_view = get_swagger_view(title='GPCRdb API')

class ProteinDetail(generics.RetrieveAPIView):
    """
    Get a single protein instance by entry name
    \n/protein/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    queryset = Protein.objects.filter(sequence_type__slug="wt").prefetch_related('family', 'species', 'source', 'residue_numbering_scheme', 'genes')
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

    queryset = ProteinFamily.objects.all().prefetch_related('parent')
    serializer_class = ProteinFamilySerializer


class ProteinFamilyDetail(generics.RetrieveAPIView):
    """
    Get a single protein family instance
    \n/proteinfamily/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    queryset = ProteinFamily.objects.all().prefetch_related("parent")
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
        queryset = ProteinFamily.objects.all().prefetch_related("parent")
        return queryset.filter(parent__slug=family)


class ProteinFamilyDescendantList(generics.ListAPIView):
    """
    Get a list of descendant families of a protein family
    \n/proteinfamily/descendants/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    serializer_class = ProteinFamilySerializer

    def get_queryset(self):
        family = self.kwargs.get('slug')
        queryset = ProteinFamily.objects.all().prefetch_related("parent")
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
        family = self.kwargs.get('slug')

        return queryset.filter(sequence_type__slug='wt', family__slug__startswith=family)\
                    .prefetch_related('family', 'species', 'source', 'residue_numbering_scheme', 'genes')


class ProteinsInFamilySpeciesList(generics.ListAPIView):
    """
    Get a list of proteins in a protein family
    \n/proteinfamily/proteins/{slug}/{species}
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{latin_name} is a species identifier from Uniprot, e.g. Homo sapiens
    """

    serializer_class = ProteinSerializer

    def get_queryset(self):
        queryset = Protein.objects.all()
        family = self.kwargs.get('slug')
        species = self.kwargs.get('latin_name')
        return queryset.filter(sequence_type__slug='wt', family__slug__startswith=family,
                               species__latin_name=species).prefetch_related('family',
                               'species', 'source', 'residue_numbering_scheme', 'genes')


class ResiduesList(generics.ListAPIView):
    """
    Get a list of residues of a protein
    \n/residues/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    serializer_class = ResidueSerializer

    def get_queryset(self):
        queryset = Residue.objects.all()
        #protein_conformation__protein__sequence_type__slug='wt',
        return queryset.filter(
            protein_conformation__protein__entry_name=self.kwargs.get('entry_name')).prefetch_related('display_generic_number','protein_segment','alternative_generic_numbers')


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
    \n{latin_name} is a species identifier from Uniprot, e.g. Homo sapiens
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

    def get(self, request, pdb_code=None, entry_name=None, representative=None):
        if pdb_code:
            structures = Structure.objects.filter(pdb_code__index=pdb_code)
        elif entry_name and representative:
            structures = Structure.objects.filter(protein_conformation__protein__parent__entry_name=entry_name,
                representative=True)
        elif entry_name:
            structures = Structure.objects.filter(protein_conformation__protein__parent__entry_name=entry_name)
        elif representative:
            structures = Structure.objects.filter(representative=True)
        else:
            structures = Structure.objects.all()

        structures = structures.exclude(refined=True).prefetch_related('protein_conformation__protein__parent__species', 'pdb_code',
            'protein_conformation__protein__parent__family', 'protein_conformation__protein__parent__species',
            'publication__web_link', 'structureligandinteraction_set__ligand__properities', 'structure_type',
            'structureligandinteraction_set__ligand__properities__ligand_type',
            'structureligandinteraction_set__ligand_role')

        # structures = self.get_structures(pdb_code, entry_name, representative)

        # convert objects to a list of dictionaries
        # normal serializers can not be used because of abstraction of tables (e.g. protein_conformation)
        s = []
        for structure in structures:
            # essential fields
            structure_data = {
                'pdb_code': structure.pdb_code.index,
                'protein': structure.protein_conformation.protein.parent.entry_name,
                'family': structure.protein_conformation.protein.parent.family.slug,
                'species': structure.protein_conformation.protein.parent.species.latin_name,
                'preferred_chain': structure.preferred_chain,
                'resolution': structure.resolution,
                'publication_date': structure.publication_date,
                'type': structure.structure_type.name,
                'state': structure.state.name,
                'distance': structure.distance,
            }

            # publication
            if structure.publication:
                structure_data['publication'] = structure.publication.web_link.__str__()
            else:
                structure_data['publication'] = None

            # ligand
            ligands = []
            for interaction in structure.structureligandinteraction_set.filter(annotated=True):
                ligand = {}
                if interaction.ligand.name:
                    ligand['name'] = interaction.ligand.name
                if interaction.ligand.properities.ligand_type and interaction.ligand.properities.ligand_type.name:
                    ligand['type'] = interaction.ligand.properities.ligand_type.name
                if interaction.ligand_role and interaction.ligand_role.name:
                    ligand['function'] = interaction.ligand_role.name
                if ligand:
                    ligands.append(ligand)
            structure_data['ligands'] = ligands

            s.append(structure_data)

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


class StructureListProtein(StructureList):
    """
    Get a list of structures of a protein
    \n/structure/protein/{entry_name}
    """


class RepresentativeStructureListProtein(StructureList):
    """
    Get a list of representative structures of a protein (one for each activation state)
    \n/structure/protein/{entry_name}/representative/
    """


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
    Get a full sequence alignment of a protein family including a consensus sequence
    \n/alignment/family/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    def get(self, request, slug=None, segments=None, latin_name=None, statistics=False):
        if slug is not None:
            # Check for specific species
            if latin_name is not None:
                ps = Protein.objects.filter(sequence_type__slug='wt', source__id=1, family__slug__startswith=slug,
                    species__latin_name=latin_name)
            else:
                ps = Protein.objects.filter(sequence_type__slug='wt', source__id=1, family__slug__startswith=slug)

            # take the numbering scheme from the first protein
            #s_slug = Protein.objects.get(entry_name=ps[0]).residue_numbering_scheme_id
            s_slug = ps[0].residue_numbering_scheme_id

            protein_family = ps[0].family.slug[:3]

            gen_list = []
            segment_list = []
            if segments is not None:
                input_list = segments.split(",")
                # fetch a list of all segments

                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True)
                for s in input_list:
                    # add to segment list
                    if s in protein_segments:
                        segment_list.append(s)
                    # get generic numbering object for generic positions
                    else:
                        # make sure the query works for all positions
                        gen_object = ResidueGenericNumberEquivalent.objects.get(label=s, scheme__id=s_slug)
                        gen_object.properties = {}
                        gen_list.append(gen_object)

                # fetch all complete protein_segments
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            if int(protein_family) < 100:
                ss = [ s for s in ss if s.proteinfamily == 'GPCR']
            elif protein_family == "100":
                ss = [ s for s in ss if s.proteinfamily == 'Gprotein']
            elif protein_family == "200":
                ss = [ s for s in ss if s.proteinfamily == 'Arrestin']

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
            if gen_list:
                a.load_segments(gen_list)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            a.calculate_statistics()

            residue_list = []
            for aa in a.full_consensus:
                residue_list.append(aa.amino_acid)

            # render the fasta template as string
            response = render_to_string('alignment/alignment_fasta.html', {'a': a}).split("\n")

            # convert the list to a dict
            ali_dict = OrderedDict({})
            for row in response:
                if row.startswith(">"):
                    k = row[1:]
                else:
                    ali_dict[k] = row
                    k = False
            ali_dict['CONSENSUS'] = ''.join(residue_list)

            # render statistics for output
            if statistics == True:
                feat = {}
                for i, feature in enumerate(AMINO_ACID_GROUPS):
                    feature_stats = a.feature_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[feature] = [item for sublist in feature_stats_clean for item in sublist]

                for i, AA in enumerate(AMINO_ACIDS):
                    feature_stats = a.amino_acid_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[AA] = [item for sublist in feature_stats_clean for item in sublist]

                ali_dict["statistics"] = feat

            return Response(ali_dict)

class FamilyAlignmentPartial(FamilyAlignment):
    """
    Get a partial sequence alignment of a protein family
    \n/alignment/family/{slug}/{segments}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    """

class FamilyAlignmentSpecies(FamilyAlignment):
    """
    Get a full sequence alignment of a protein family
    \n/alignment/family/{slug}//{species}
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{species} is a species identifier from Uniprot, e.g. Homo sapiens
    """

class FamilyAlignmentPartialSpecies(FamilyAlignment):
    """
    Get a partial sequence alignment of a protein family
    \n/alignment/family/{slug}/{segments}/{species}
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    \n{species} is a species identifier from Uniprot, e.g. Homo sapiens
    """


class ProteinSimilaritySearchAlignment(views.APIView):
    """
    Get a segment sequence alignment of two or more proteins ranked by similarity
    \n/alignment/similarity/{proteins}/{segments}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human,cxcr4_human,
    where the first protein is the query protein and the following the proteins to compare it to
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    """

    def get(self, request, proteins=None, segments=None):
        if proteins is not None:
            protein_list = proteins.split(",")
            # first in API should be reference
            ps = Protein.objects.filter(sequence_type__slug='wt', entry_name__in=protein_list[1:])
            reference = Protein.objects.filter(sequence_type__slug='wt', entry_name__in=[protein_list[0]])

            # take the numbering scheme from the first protein
            s_slug = Protein.objects.get(entry_name=protein_list[0]).residue_numbering_scheme_id

            protein_family = ps[0].family.slug[:3]

            gen_list = []
            segment_list = []
            if segments is not None:
                input_list = segments.split(",")
                # fetch a list of all segments
                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True)
                for s in input_list:
                    # add to segment list
                    if s in protein_segments:
                        segment_list.append(s)
                    # get generic numbering object for generic positions
                    else:
                        # make sure the query works for all positions
                        gen_object = ResidueGenericNumberEquivalent.objects.get(label=s, scheme__id=s_slug)
                        gen_object.properties = {}
                        gen_list.append(gen_object)

                # fetch all complete protein_segments
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            if int(protein_family) < 100:
                ss = [ s for s in ss if s.proteinfamily == 'GPCR']
            elif protein_family == "100":
                ss = [ s for s in ss if s.proteinfamily == 'Gprotein']
            elif protein_family == "200":
                ss = [ s for s in ss if s.proteinfamily == 'Arrestin']

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from API into the alignment
            a.load_reference_protein(reference[0])
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
            if gen_list:
                a.load_segments(gen_list)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            # calculate identity and similarity of each row compared to the reference
            a.calculate_similarity()

            # render the fasta template as string
            response = render_to_string('alignment/alignment_fasta.html', {'a': a}).split("\n")

            # convert the list to a dict
            ali_dict = {}
            k = False
            num = 0
            for i, row in enumerate(response):
                if row.startswith(">"):
                    k = row[1:]
                elif k:
                    # add the query as 100 identical/similar to the beginning (like on the website)
                    if num == 0:
                        a.proteins[num].identity = 100
                        a.proteins[num].similarity = 100
                    # order dict after custom list
                    keyorder = ["similarity","identity","AA"]
                    ali_dict[k] = {"AA": row, "identity": int(str(a.proteins[num].identity).replace(" ","")),
                    "similarity": int(str(a.proteins[num].similarity).replace(" ",""))}
                    ali_dict[k] = OrderedDict(sorted(ali_dict[k].items(), key=lambda t: keyorder.index(t[0])))
                    num+=1
                    k = False
            ali_dict_ordered = OrderedDict(sorted(ali_dict.items(), key=lambda x: x[1]['similarity'], reverse=True))
            return Response(ali_dict_ordered)

class ProteinAlignment(views.APIView):
    """
    Get a full sequence alignment of two or more proteins
    \n/alignment/protein/{proteins}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human
    """

    def get(self, request, proteins=None, segments=None, statistics=False):
        if proteins is not None:
            protein_list = proteins.split(",")
            ps = Protein.objects.filter(sequence_type__slug='wt', entry_name__in=protein_list)

            # take the numbering scheme from the first protein
            #s_slug = Protein.objects.get(entry_name=protein_list[0]).residue_numbering_scheme_id
            s_slug = ps[0].residue_numbering_scheme_id

            protein_family = ps[0].family.slug[:3]

            gen_list = []
            segment_list = []
            if segments is not None:
                input_list = segments.split(",")

                # fetch a list of all segments
                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True)
                for s in input_list:
                    # add to segment list
                    if s in protein_segments:
                        segment_list.append(s)
                    # get generic numbering object for generic positions
                    else:
                        gen_object = ResidueGenericNumberEquivalent.objects.get(label=s, scheme__id=s_slug)
                        gen_object.properties = {}
                        gen_list.append(gen_object)

                # fetch all complete protein_segments
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            if int(protein_family) < 100:
                ss = [ s for s in ss if s.proteinfamily == 'GPCR']
            elif protein_family == "100":
                ss = [ s for s in ss if s.proteinfamily == 'Gprotein']
            elif protein_family == "200":
                ss = [ s for s in ss if s.proteinfamily == 'Arrestin']

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
            if gen_list:
                a.load_segments(gen_list)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            # calculate statistics
            if statistics == True:
                a.calculate_statistics()

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

            # render statistics for output
            if statistics == True:
                feat = {}
                for i, feature in enumerate(AMINO_ACID_GROUPS):
                    feature_stats = a.feature_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[feature] = [item for sublist in feature_stats_clean for item in sublist]

                for i, AA in enumerate(AMINO_ACIDS):
                    feature_stats = a.amino_acid_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[AA] = [item for sublist in feature_stats_clean for item in sublist]

                ali_dict["statistics"] = feat

            return Response(ali_dict)

class ProteinAlignmentStatistics(ProteinAlignment):
    """
    Add a /statics at the end of an alignment in order to
    receive an additional residue property statistics output e.g.:
    \n/alignment/protein/{proteins}/{segments}/statistics
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    """

class ProteinAlignmentPartial(ProteinAlignment):
    """
    Get a partial sequence alignment of two or more proteins
    \n/alignment/protein/{proteins}/{segments}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
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
                input_list = segments.split(",")
                ss = ProteinSegment.objects.filter(slug__in=input_list, partial=False)
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


class StructureAssignGenericNumbers(views.APIView):
    """
    Assign generic residue numbers (Ballesteros-Weinstein and GPCRdb schemes) to an uploaded pdb file.
    \n/structure/assign_generic_numbers\n
    e.g.
    curl -X POST -F "pdb_file=@myfile.pdb" http://gpcrdb.org/services/structure/assign_generic_numbers
    """
    parser_classes = (FileUploadParser,)
    renderer_classes = (PDBRenderer, )

    def post(self, request):

        # root, ext = os.path.splitext(request._request.FILES['pdb_file'].name)
        generic_numbering = GenericNumbering(StringIO(request._request.FILES['pdb_file'].file.read().decode('UTF-8',"ignore")))
        out_struct = generic_numbering.assign_generic_numbers()
        out_stream = StringIO()
        io = PDBIO()
        io.set_structure(out_struct)
        io.save(out_stream)
        print(len(out_stream.getvalue()))
        # filename="{}_GPCRdb.pdb".format(root)
        return Response(out_stream.getvalue())


class StructureSequenceParser(views.APIView):
    """
    Analyze the uploaded pdb structure listing auxiliary proteins, mutations, deletions and insertions.
    \n/structure/structure/parse_pdb\n
    e.g.
    curl -X POST -F "pdb_file=@myfile.pdb" http://gpcrdb.org/services/structure/parse_pdb
    """
    parser_classes = (FileUploadParser,)
    renderer_classes = (JSONRenderer, )

    def post(self, request):
        # root, ext = os.path.splitext(request._request.FILES['pdb_file'].name)
        pdb_file = StringIO(request._request.FILES['pdb_file'].file.read().decode('UTF-8',"ignore"))
        header = parse_pdb_header(pdb_file)
        parser = SequenceParser(pdb_file)

        json_data = OrderedDict()
        json_data["header"] = header
        json_data.update(parser.get_fusions())
        json_data.update(parser.get_mutations())
        json_data.update(parser.get_deletions())

        return Response(json_data)


class StructureLigandInteractions(generics.ListAPIView):
    """
    Get a list of interactions between structure and ligand
    \n/structure/{pdb_code}/interaction/
    \n{pdb_code} is a structure identifier from the Protein Data Bank, e.g. 2RH1
    """
    serializer_class = StructureLigandInteractionSerializer

    def get_queryset(self):
        queryset = ResidueFragmentInteraction.objects.all()
        queryset = queryset.prefetch_related('structure_ligand_pair__structure__pdb_code',
                                             'interaction_type',
                                             'fragment__residue__generic_number',
                                             'fragment__residue__display_generic_number',
                                             )
        #queryset = queryset.exclude(interaction_type__type='hidden').order_by('fragment__residue__sequence_number')
        queryset = queryset.order_by('fragment__residue__sequence_number')
        slug = self.kwargs.get('pdb_code')
        return queryset.filter(structure_ligand_pair__structure__pdb_code__index=slug,
                               structure_ligand_pair__annotated=True)


class MutantList(generics.ListAPIView):
    """
    Get a list of mutants of single protein instance by entry name
    \n/mutant/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """
    serializer_class = MutationSerializer

    def get_queryset(self):
        queryset = MutationRaw.objects.all()
        return queryset.filter(protein=self.kwargs.get('entry_name'))

class DrugList(views.APIView):
    """
    Get a list of drugs for a single protein instance by entry name
    \n/drugs/{proteins}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """
    def get(self, request, entry_name=None):

        drugs = Drugs.objects.filter(target__entry_name=entry_name).distinct()

        druglist = []
        for drug in drugs:
            drugname = drug.name
            drugtype = drug.drugtype
            clinical = drug.clinicalstatus
            phasedate = drug.phasedate
            if clinical != '-':
                status = drug.status + ' (' + drug.clinicalstatus + ', ' + phasedate + ')'
            else:
                status = drug.status
            approval = drug.approval
            indication = drug.indication
            moa = drug.moa
            novelty = drug.novelty
            druglist.append({'name':drugname, 'approval': approval, 'indication': indication, 'status':status, 'drugtype':drugtype, 'moa':moa, 'novelty': novelty})

        return Response(druglist)
