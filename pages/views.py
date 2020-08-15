from django.shortcuts import render

from common.models import ReleaseNotes


def releasenotes(request):
    """Get release notes"""
    context = {}
    context['release_notes'] = ReleaseNotes.objects.all()
    return render(request, 'pages/releasenotes.html', context)
