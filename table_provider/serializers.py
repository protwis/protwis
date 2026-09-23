from rest_framework import serializers
from common.models import WebResource, WebLink

class GpcrStatisticsSummarySerializer(serializers.Serializer):
    """Serializer to generate JSON for the structure statistics data for GPCRs stored in the gpcr_structure_statistics database table."""

    _uniprot_resource = None
    _entrez_resource = None

    uniprot = serializers.SerializerMethodField()
    protein_name = serializers.CharField(read_only=True)
    receptor_family = serializers.ReadOnlyField()
    receptor_class = serializers.ReadOnlyField(read_only=True, )
    gene = serializers.SerializerMethodField()
    species_common_name = serializers.CharField(read_only=True)

    # Structure counts by state
    st_inact_c = serializers.ReadOnlyField(read_only=True)
    st_interm_c = serializers.ReadOnlyField(read_only=True)
    st_act_c = serializers.ReadOnlyField(read_only=True)

    # Best resolution by state ("-" when no structure of that state exists)
    best_res_inact = serializers.SerializerMethodField()
    best_res_interm = serializers.SerializerMethodField()
    best_res_act = serializers.SerializerMethodField()

    # Transducer / extra-protein counts (map to annotated names)
    transd_gio_c = serializers.ReadOnlyField(read_only=True)
    transd_gq11_c = serializers.ReadOnlyField(read_only=True)
    transd_g1213_c = serializers.ReadOnlyField(read_only=True)
    transd_arrest_c = serializers.ReadOnlyField(read_only=True)
    transd_grks_c = serializers.ReadOnlyField(read_only=True)

    # Ligand type counts
    ligand_types_smmol_c = serializers.ReadOnlyField(read_only=True)
    ligand_types_pep_c = serializers.ReadOnlyField(read_only=True)

    # Modality group counts (map to annotated names)
    limodgrp_stim_c = serializers.ReadOnlyField(read_only=True)
    limodgrp_inhib_c = serializers.ReadOnlyField(read_only=True)
    limodgrp_unk_c = serializers.ReadOnlyField(read_only=True)

    # Individual modality counts
    limod_apo_c = serializers.ReadOnlyField(read_only=True)
    limod_ag_c = serializers.ReadOnlyField(read_only=True)
    limod_ag_partial_c = serializers.ReadOnlyField(read_only=True)
    limod_allo_antag_c = serializers.ReadOnlyField(read_only=True)
    limod_antag_c = serializers.ReadOnlyField(read_only=True)
    limod_inv_ag_c = serializers.ReadOnlyField(read_only=True)
    limod_nam_c = serializers.ReadOnlyField(read_only=True)
    limod_pam_c = serializers.ReadOnlyField(read_only=True)
    limod_unk_c = serializers.ReadOnlyField(read_only=True)

    class Meta:
        """The Meta data subclass provides additional information used for filter processing.

        Meta.datatypes: contains data types of each field which is used when  building server-side filters
        Meta.serializer_method_to_filter_field_map: contains a mapping for fields that use SerializerMethodField (e.g. weblinks) to a single source
        field in the database (serializer_method_to_filter_field_map).
        """

        datatypes = {
            'uniprot' : 'text', 'protein_name' : 'text', 'receptor_family' : 'text', 'receptor_class' : 'text', 'gene' : 'text', 'species_common_name' : 'text',
            'st_inact_c' : 'numeric', 'st_interm_c' : 'numeric', 'st_act_c' : 'numeric',
            'best_res_inact' : 'numeric', 'best_res_interm' : 'numeric', 'best_res_act' : 'numeric',
            'transd_gio_c': 'numeric', 'transd_gq11_c': 'numeric', 'transd_g1213_c': 'numeric', 'transd_arrest_c': 'numeric', 'transd_grks_c': 'numeric',
            'ligand_types_smmol_c': 'numeric', 'ligand_types_pep_c': 'numeric',
            'limodgrp_stim_c': 'numeric', 'limodgrp_inhib_c': 'numeric', 'limodgrp_unk_c': 'numeric', 'limod_apo_c': 'numeric', 'limod_ag_c': 'numeric',
            'limod_ag_partial_c': 'numeric', 'limod_allo_antag_c': 'numeric', 'limod_antag_c': 'numeric', 'limod_inv_ag_c': 'numeric', 'limod_nam_c': 'numeric',
            'limod_pam_c': 'numeric', 'limod_unk_c': 'numeric',
        }
        #Mapping for fields that use a SerializerMethodField to a single source field in the database. This is used for filtering and sorting.
        serializer_method_to_filter_field_map = {
            'uniprot': 'uniprot_entry_name', 'gene': 'gene_name',
            'best_res_inact': 'best_res_inact', 'best_res_interm': 'best_res_interm', 'best_res_act': 'best_res_act',
        }

    @staticmethod
    def _format_resolution(value):
        return value if value is not None else "-"

    def get_best_res_inact(self, obj):
        return self._format_resolution(getattr(obj, 'best_res_inact', None))

    def get_best_res_interm(self, obj):
        return self._format_resolution(getattr(obj, 'best_res_interm', None))

    def get_best_res_act(self, obj):
        return self._format_resolution(getattr(obj, 'best_res_act', None))

    def get_gene_weblink(self, entrez_id):
        if entrez_id:
            return str(WebLink(index=entrez_id, web_resource=self._get_entrez_resource()))
        return None

    def get_gene(self, obj):
        if hasattr(obj, 'gene_name') and obj.gene_name:
            gene_anchor = {}
            gene_anchor['anchor_content'] = obj.gene_name
            if hasattr(obj, 'gene_entrez_id') and obj.gene_entrez_id:
                gene_anchor['anchor_href'] = self.get_gene_weblink(obj.gene_entrez_id)
            else:
                gene_anchor['anchor_href'] = None
            return gene_anchor
        return None

    def get_uniprot_weblink(self, obj):
        if obj.uniprot_accession:
            return str(WebLink(index=obj.uniprot_accession, web_resource=self._get_uniprot_resource()))
        return None

    def get_uniprot(self, obj):
        uniprot_anchor = {}
        uniprot_anchor['anchor_content'] = obj.uniprot_entry_name
        uniprot_anchor['anchor_href'] = self.get_uniprot_weblink(obj)
        return uniprot_anchor

    #Make WebResource objects for Uniprot and Entrez Gene available at the class level to prevent multiple database queries.
    def _get_uniprot_resource(self):
        if self._uniprot_resource is None:
            self._uniprot_resource = WebResource.objects.get(slug="uniprot")
        return self._uniprot_resource

    def _get_entrez_resource(self):
        if self._entrez_resource is None:
            self._entrez_resource = WebResource.objects.get(slug="entrez_gene")
        return self._entrez_resource

    @classmethod
    def source_mapping(cls):
        """Generate a mapping of serializer field names to their corresponding source fields in the database.

        Handles both simple fields and those that use SerializerMethodField, ensuring that all fields are accounted for in the mapping.
            :return: A dictionary mapping serializer field names to database source field names.
        """
        mapping = dict()
        for key, field in cls.__dict__['_declared_fields'].items():

            if isinstance(field, serializers.SerializerMethodField):
                try:
                    mapping[key] = cls.Meta.serializer_method_to_filter_field_map[key]
                    continue
                except KeyError:
                    raise KeyError(f"No query source mapping found for field '{key}' in serializer. " + \
                                    "Ensure you have declared all fields that do not have a source " + \
                                    "parameter in serializer_method_to_filter_field_map under meta data")

            mapping[key] = key

        return mapping


