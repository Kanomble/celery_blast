from django.http import Http404
from django.shortcuts import redirect

#redirects to the login page if user is not logged in
def unauthenticated_user(view_func):
    def wrapper_func(request,*args,**kwargs):
        if request.user.is_authenticated:
            return redirect('main')
        else:
            return view_func(request,*args,**kwargs)
    return wrapper_func

#TODO use this decorator at appropriate functions
def allowed_user(allowed_roles=[]):
    def decorator(view_func):
        def wrapper_func(request,*args,**kwargs):
            #print("Working",allowed_roles)
            group = None
            if request.user.groups.exists():
                group = request.user.groups.all()[0].name
            if group in allowed_roles:
                return view_func(request,*args,**kwargs)
            else:
                return Http404
        return wrapper_func
    return decorator