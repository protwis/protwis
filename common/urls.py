from django.conf.urls import url, include

from common import views


urlpatterns = [
    url(r'^addtoselection', views.AddToSelection, name='addtoselection'),
    url(r'^removefromselection', views.RemoveFromSelection, name='removefromselection'),
    url(r'^clearselection', views.ClearSelection, name='clearselection'),
    url(r'^selectrange', views.SelectRange, name='selectrange'),
    url(r'^togglefamilytreenode', views.ToggleFamilyTreeNode, name='togglefamilytreenode'),
    url(r'^selectionannotation', views.SelectionAnnotation, name='selectionannotation'),
    url(r'^selectionspeciespredefined', views.SelectionSpeciesPredefined, name='selectionspeciespredefined'),
    url(r'^selectionspeciestoggle', views.SelectionSpeciesToggle, name='selectionspeciestoggle'),
    url(r'^expandsegment', views.ExpandSegment, name='expandsegment'),
    url(r'^selectfullsequence', views.SelectFullSequence, name='selectfullsequence'),
    url(r'^selectalignablesegments', views.SelectAlignableSegments, name='selectalignablesegments'),
    url(r'^selectalignableresidues', views.SelectAlignableResidues, name='selectalignableresidues'),
    url(r'selectionschemespredefined', views.SelectionSchemesPredefined, name='selectionschemespredefined'),
    url(r'selectionschemestoggle', views.SelectionSchemesToggle, name='selectionschemestoggle'),
    url(r'settreeselection', views.SetTreeSelection, name='settreeselection'),
    url(r'selectresiduefeature', views.SelectResidueFeature, name='selectresiduefeature'),
    url(r'addresiduegroup', views.AddResidueGroup, name='addresiduegroup'),
    url(r'selectresiduegroup', views.SelectResidueGroup, name='selectresiduegroup'),
    url(r'removeresiduegroup', views.RemoveResidueGroup, name='removeresiduegroup'),
    url(r'setgroupminmatch', views.SetGroupMinMatch, name='setgroupminmatch'),
    url(r'residuesdownload', views.ResiduesDownload, name='residuesupload'),
    url(r'residuesupload', views.ResiduesUpload, name='residuesupload'),
    url(r'^selectiongproteinpredefined', views.SelectionGproteinPredefined, name='selectiongproteinpredefined'),
    url(r'^selectiongproteintoggle', views.SelectionGproteinToggle, name='selectiongproteintoggle'),
    url(r'^targetformread', views.ReadTargetInput, name='targetformread'),
    url(r'^exportexcelsuggestions$', views.ExportExcelSuggestions, name='exportexcelsuggestions'),
    url(r'^exportexcelmodifications$', views.ExportExcelModifications, name='exportexcelmodifications'),
    url(r'^exportexceldownload/(?P<ts>[^/]*?)/(?P<entry_name>.+)$', views.ExportExcelDownload, name='exportexceldownload'),
    url(r'^importexcel$', views.ImportExcel, name='importexcel'),
]