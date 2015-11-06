from common.models import ReleaseNotes

def releasenotes(request):
    context = {}
    context['release_notes'] = ReleaseNotes.objects.all()
    return render(request, 'pages/releasenotes.html', context)
