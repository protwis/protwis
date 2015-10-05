from django.conf.urls import patterns, url
from structure.views import *
from structure import views
from django.conf import settings
from django.views.generic import TemplateView

urlpatterns = patterns('',
    url(r'^$', StructureBrowser.as_view(), name='structure_browser'),
    url(r'^selection_convert$', ConvertStructuresToProteins, name='convert'),
    url(r'^template_browser', TemplateBrowser.as_view(), name='structure_browser'),
    url(r'^template_selection', TemplateTargetSelection.as_view(), name='structure_browser'),
    url(r'^template_segment_selection', TemplateSegmentSelection.as_view(), name='structure_browser'),
    url(r'^statistics$', StructureStatistics.as_view(), name='structure_statistics'),
    url(r'pdb_download', PDBDownload.as_view(), name='pdb_download'),
    url(r'^generic_numbering_index', GenericNumberingIndex.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_results$', GenericNumberingResults.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_results/(?P<substructure>\w+)$', GenericNumberingDownload.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_selection', GenericNumberingSelection.as_view(), name='generic_numbering'),
    url(r'^superposition_workflow_index$', SuperpositionWorkflowIndex.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_index/(?P<clear>\w{4})$', SuperpositionWorkflowIndex.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_selection', SuperpositionWorkflowSelection.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_results$', SuperpositionWorkflowResults.as_view(), name='superposition_workflow'),
    url(r'^superposition_workflow_results/(?P<substructure>\w+)$', SuperpositionWorkflowDownload.as_view(), name='superposition_workflow'),
    url(r'^fragment_superposition_index', FragmentSuperpositionIndex.as_view(), name='fragment_superposition'),
    url(r'^fragment_superposition_results', FragmentSuperpositionResults.as_view(), name='fragment_superposition'),
    url(r'^output/(?P<outfile>\w+.\w{3})/(?P<replacement_tag>\w+)$', ServePdbOutfile, name='structural_tools_result'),
    url(r'^zipoutput/(?P<outfile>\w+.\w{3})/', ServeZipOutfile, name='structural_tools_result'), 
    url(r'^showtrees', RenderTrees, name='render'),
    url(r'^(?P<pdbname>\w+)$', StructureDetails, name='structure_details'), 
    url(r'^pdb/(?P<pdbname>\w+)$', ServePdbDiagram, name='structure_serve_pdb'), 
    url(r'^pdb/(?P<pdbname>\w+)/ligand/(?P<ligand>.+)$', ServePdbLigandDiagram, name='structure_serve_pdb_ligand'), 

)