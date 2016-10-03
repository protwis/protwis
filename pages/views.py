from common.models import ReleaseNotes
from django.shortcuts import get_object_or_404, render

def releasenotes(request):
    context = {}
    context['release_notes'] = ReleaseNotes.objects.all()
    return render(request, 'pages/releasenotes.html', context)
