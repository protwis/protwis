from django.conf.urls import url

from seqsign import views

urlpatterns = [
    url(r'^$', (views.PosTargetSelection.as_view()), name='pgselection'),
    url(r'^segmentselectionsignature', views.SegmentSelectionSignature.as_view(), name='segmentselectionsignature'),
    url(r'^render_signature_match_scores/(?P<cutoff>[\d]+)', views.render_signature_match_scores, name='render_signature_match_scores'),
    url(r'render_signature_match_excel', views.render_signature_match_excel, name=''),
    url(r'^render_signature_excel', views.render_signature_excel, name='render_signature_excel'),
    url(r'^render_signature', views.render_signature, name='rendersignature'),
    url(r'render_positive', views.render_reordered, {'group' : 'positive'}, name='render-reordered'),
    url(r'render_negative', views.render_reordered, {'group' : 'negative'}, name='render-reordered'),
    url(r'savepos', views.preserve_targets, name='ngselection'),
    url(r'negativegroupselection', views.NegTargetSelection.as_view(), name='ngselection'),
]
