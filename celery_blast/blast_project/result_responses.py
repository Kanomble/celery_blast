from django.http import HttpResponse
from django.utils.cache import patch_vary_headers


HTML_CONTENT_TYPE = "text/html; charset=utf-8"

STATIC_RESULT_CSP = (
    "default-src 'none'; "
    "img-src 'self' data:; "
    "style-src 'self' 'unsafe-inline'; "
    "font-src 'self' data:; "
    "object-src 'none'; "
    "base-uri 'none'; "
    "form-action 'none'; "
    "frame-ancestors 'none'; "
    "sandbox allow-downloads"
)

INTERACTIVE_RESULT_CSP = (
    "default-src 'none'; "
    "script-src 'self' 'unsafe-inline' 'unsafe-eval' "
    "https://code.jquery.com https://cdn.datatables.net https://cdn.bokeh.org "
    "https://cdn.jsdelivr.net https://cdnjs.cloudflare.com; "
    "style-src 'self' 'unsafe-inline' "
    "https://cdn.datatables.net https://cdn.bokeh.org https://cdn.jsdelivr.net "
    "https://cdnjs.cloudflare.com; "
    "img-src 'self' data: blob:; "
    "font-src 'self' data: https://cdnjs.cloudflare.com; "
    "connect-src 'self' data: blob:; "
    "media-src 'self' data: blob:; "
    "object-src 'none'; "
    "base-uri 'none'; "
    "form-action 'none'; "
    "frame-ancestors 'none'; "
    "sandbox allow-scripts allow-forms allow-popups allow-downloads"
)


def _coerce_html(html_data):
    if isinstance(html_data, (list, tuple)):
        return "".join(str(item) for item in html_data)
    return str(html_data)


def is_html_content_type(content_type):
    return content_type.split(";", 1)[0].strip().lower() in {
        "text/html",
        "application/xhtml+xml",
    }


def apply_result_security_headers(response, interactive=False):
    response["X-Content-Type-Options"] = "nosniff"
    response["X-Frame-Options"] = "DENY"
    response["Content-Security-Policy"] = INTERACTIVE_RESULT_CSP if interactive else STATIC_RESULT_CSP
    patch_vary_headers(response, ("Cookie",))
    return response


def generated_html_response(html_data, interactive=False):
    response = HttpResponse(_coerce_html(html_data), content_type=HTML_CONTENT_TYPE)
    return apply_result_security_headers(response, interactive=interactive)
