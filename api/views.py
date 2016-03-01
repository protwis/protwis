from django.shortcuts import render
from rest_framework import views, generics, viewsets
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser, FormParser, FileUploadParser
from django.template.loader import render_to_string
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinFamily, Species, ProteinSegment
from residue.models import Residue, ResidueGenericNumber, ResidueNumberingScheme, ResidueGenericNumberEquivalent
from structure.models import Structure
from structure.assign_generic_numbers_gpcr import GenericNumbering
from api.serializers import (ProteinSerializer, ProteinFamilySerializer, SpeciesSerializer, ResidueSerializer,
    ResidueExtendedSerializer, StructureSerializer)
from api.renderers import PDBRenderer
from common.alignment import Alignment

import json, os
from io import StringIO
from Bio.PDB import PDBIO
from collections import OrderedDict

# FIXME add
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
        return queryset.filter(sequence_type__slug='wt', family__slug__startswith=self.kwargs.get('slug'), 
            species__latin_name="Homo sapiens")


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

        structures = structures.prefetch_related('protein_conformation__protein__parent__species', 'pdb_code',
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

    def get(self, request, slug=None, segments=None, latin_name=None):

        if slug is not None:
            # Check for specific species
            if latin_name is not None:
                ps = Protein.objects.filter(sequence_type__slug='wt', source__id=1, family__slug__startswith=slug, 
                    species__latin_name=latin_name)
            else:
                ps = Protein.objects.filter(sequence_type__slug='wt', source__id=1, family__slug__startswith=slug)
            
            # take the numbering scheme from the first protein
            s_slug = Protein.objects.get(entry_name=ps[0]).residue_numbering_scheme_id

            if segments is not None:
                input_list = segments.split(",")
                # fetch a list of all segments
                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True) 
                gen_list = []
                segment_list = []
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

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
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
            return Response(ali_dict)


class FamilyAlignmentPartial(FamilyAlignment):
    """
    Get a partial sequence alignment of a protein family
    \n/alignment/family/{slug}/{segments}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{segments} is a comma separated list of protein segment identifiers and/ or 
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
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
    \n/alignment/similarity/{proteins}/
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

            if segments is not None:
                input_list = segments.split(",")
                # fetch a list of all segments
                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True) 
                gen_list = []
                segment_list = []
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

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from API into the alignment
            a.load_reference_protein(reference)
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
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

    def get(self, request, proteins=None, segments=None):
        if proteins is not None:
            protein_list = proteins.split(",")
            ps = Protein.objects.filter(sequence_type__slug='wt', entry_name__in=protein_list)

            # take the numbering scheme from the first protein
            s_slug = Protein.objects.get(entry_name=protein_list[0]).residue_numbering_scheme_id

            if segments is not None:
                input_list = segments.split(",")
                # fetch a list of all segments
                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True) 
                gen_list = []
                segment_list = []
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

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
            a.load_segments(gen_list)
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

        root, ext = os.path.splitext(request.FILES['pdb_file'].name)
        generic_numbering = GenericNumbering(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8',"ignore")))
        out_struct = generic_numbering.assign_generic_numbers()
        out_stream = StringIO()
        io = PDBIO()
        io.set_structure(out_struct)
        io.save(out_stream)
        print(len(out_stream.getvalue()))
        # filename="{}_GPCRdb.pdb".format(root)
        return Response(out_stream.getvalue())
