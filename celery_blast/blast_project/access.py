from functools import wraps

from django.shortcuts import get_object_or_404

from .models import BlastProject, RemoteBlastProject


LOCAL_PROJECT = 'local'
REMOTE_PROJECT = 'remote'


def normalize_project_scope(remote_or_local):
    if remote_or_local == LOCAL_PROJECT:
        return LOCAL_PROJECT
    if remote_or_local == REMOTE_PROJECT:
        return REMOTE_PROJECT
    raise ValueError("project scope must be 'local' or 'remote'")


def get_owned_project_or_404(user, project_id, remote_or_local=LOCAL_PROJECT):
    scope = normalize_project_scope(remote_or_local)
    if scope == LOCAL_PROJECT:
        return get_object_or_404(BlastProject, id=project_id, project_user=user)
    return get_object_or_404(RemoteBlastProject, id=project_id, r_project_user=user)


def get_owned_local_project_or_404(user, project_id):
    return get_owned_project_or_404(user, project_id, LOCAL_PROJECT)


def get_owned_remote_project_or_404(user, project_id):
    return get_owned_project_or_404(user, project_id, REMOTE_PROJECT)


def project_owner_required(remote_or_local=None, project_id_kwarg='project_id'):
    def decorator(view_func):
        @wraps(view_func)
        def wrapped(request, *args, **kwargs):
            project_id = kwargs[project_id_kwarg]
            scope = remote_or_local
            if scope is None:
                scope = kwargs.get('remote_or_local', LOCAL_PROJECT)
            request.owned_blast_project = get_owned_project_or_404(request.user, project_id, scope)
            return view_func(request, *args, **kwargs)

        return wrapped

    return decorator
