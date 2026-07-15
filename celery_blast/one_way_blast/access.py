from functools import wraps

from django.shortcuts import get_object_or_404

from .models import OneWayBlastProject, OneWayRemoteBlastProject


def get_owned_one_way_project_or_404(user, project_id, remote=False):
    if remote:
        return get_object_or_404(OneWayRemoteBlastProject, id=project_id, r_project_user=user)
    return get_object_or_404(OneWayBlastProject, id=project_id, project_user=user)


def one_way_project_owner_required(remote=None):
    def decorator(view_func):
        @wraps(view_func)
        def wrapped(request, *args, **kwargs):
            is_remote = remote
            if is_remote is None:
                if 'remote' in kwargs:
                    is_remote = bool(kwargs['remote'])
                else:
                    is_remote = kwargs.get('project_type') == 'remote_searches'
            request.owned_one_way_project = get_owned_one_way_project_or_404(
                request.user,
                kwargs['project_id'],
                is_remote,
            )
            return view_func(request, *args, **kwargs)

        return wrapped

    return decorator
