import mimetypes
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
from urllib.parse import quote, unquote

from django.conf import settings
from django.http import FileResponse, Http404, HttpResponse
from django.utils.http import content_disposition_header

from .result_responses import apply_result_security_headers, is_html_content_type


def _apply_protected_file_headers(response, is_html):
    if is_html:
        apply_result_security_headers(response, interactive=True)
    else:
        response['X-Content-Type-Options'] = 'nosniff'
    return response


def absolute_project_path(path):
    path = Path(path)
    if path.is_absolute():
        return path
    return Path(settings.BASE_DIR) / path


def _decode_relative_path(relative_path):
    decoded = str(relative_path)
    for _ in range(3):
        next_value = unquote(decoded)
        if next_value == decoded:
            break
        decoded = next_value
    return decoded


def _validate_relative_path(relative_path):
    decoded = _decode_relative_path(relative_path)
    posix_path = PurePosixPath(decoded)
    windows_path = PureWindowsPath(decoded)
    if (
        not decoded
        or '\x00' in decoded
        or '\\' in decoded
        or ':' in decoded
        or posix_path.is_absolute()
        or windows_path.is_absolute()
        or any(part in ('', '.', '..') for part in posix_path.parts)
    ):
        raise Http404("File not found")
    return posix_path


def resolve_project_file(project_root, relative_path):
    safe_path = _validate_relative_path(relative_path)
    try:
        root = absolute_project_path(project_root).resolve(strict=True)
    except OSError:
        raise Http404("File not found")

    candidate = (root / Path(*safe_path.parts)).resolve(strict=False)
    try:
        if os.path.commonpath([str(root), str(candidate)]) != str(root):
            raise Http404("File not found")
    except ValueError:
        raise Http404("File not found")

    if not candidate.is_file():
        raise Http404("File not found")
    return candidate


def _x_accel_redirect_path(file_path):
    try:
        media_root = absolute_project_path(settings.PROTECTED_MEDIA_ROOT).resolve(strict=True)
        relative_path = file_path.resolve(strict=True).relative_to(media_root)
    except (OSError, ValueError):
        return None
    internal_url = settings.PROTECTED_MEDIA_INTERNAL_URL.rstrip('/') + '/'
    return internal_url + quote(relative_path.as_posix(), safe='/')


def protected_file_response(file_path, as_attachment=False):
    content_type, encoding = mimetypes.guess_type(str(file_path))
    content_type = content_type or 'application/octet-stream'
    filename = file_path.name
    is_html = is_html_content_type(content_type)

    use_x_accel = (
        settings.PROTECTED_MEDIA_USE_X_ACCEL
        and file_path.stat().st_size >= settings.PROTECTED_MEDIA_X_ACCEL_SIZE_THRESHOLD
    )
    if use_x_accel:
        redirect_path = _x_accel_redirect_path(file_path)
        if redirect_path is not None:
            response = HttpResponse(content_type=content_type)
            response['X-Accel-Redirect'] = redirect_path
            response['Content-Disposition'] = content_disposition_header(as_attachment, filename)
            if encoding:
                response['Content-Encoding'] = encoding
            _apply_protected_file_headers(response, is_html)
            return response

    response = FileResponse(
        open(file_path, 'rb'),
        as_attachment=as_attachment,
        filename=filename,
        content_type=content_type,
    )
    _apply_protected_file_headers(response, is_html)
    return response


def serve_protected_project_file(project_root, relative_path, as_attachment=False):
    file_path = resolve_project_file(project_root, relative_path)
    return protected_file_response(file_path, as_attachment=as_attachment)
