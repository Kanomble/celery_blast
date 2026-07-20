"""Static checks for Docker build policy."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCKERFILES = [ROOT / "Dockerfile", ROOT / "nginx" / "Dockerfile"]


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def fail(message: str) -> None:
    raise AssertionError(message)


def dockerfile_text() -> str:
    return "\n".join(read(path) for path in DOCKERFILES)


def assert_base_images_are_explicit() -> None:
    for path in DOCKERFILES:
        for line in read(path).splitlines():
            stripped = line.strip()
            if not stripped.startswith("FROM "):
                continue
            if "${UBUNTU_FOCAL_AMD64_IMAGE}" in stripped:
                continue
            if stripped.endswith(":latest") or stripped == "FROM nginx":
                fail(f"{path} must use an explicit non-latest base image tag: {stripped}")

    dockerfile = read(ROOT / "Dockerfile")
    if "UBUNTU_FOCAL_AMD64_IMAGE=ubuntu:focal" not in dockerfile:
        fail("Dockerfile must declare ubuntu:focal through UBUNTU_FOCAL_AMD64_IMAGE")


def assert_no_remote_shell_execution() -> None:
    text = dockerfile_text()
    forbidden_patterns = [
        r"curl\s+[^|\n]+\|\s*(sh|bash|tar)",
        r"wget\s+[^|\n]+\|\s*(sh|bash|tar)",
        r"sh\s+-c\s+\"\$\(curl",
        r"bash\s+-c\s+\"\$\(curl",
        r"ftp://",
        r"refs/heads/trimAl",
        r"conda update --all",
    ]
    for pattern in forbidden_patterns:
        if re.search(pattern, text):
            fail(f"Dockerfiles contain forbidden remote acquisition pattern: {pattern}")

    if "--proto '=https' --tlsv1.2" not in text:
        fail("Dockerfile downloads must constrain curl to HTTPS")


def assert_remote_downloads_are_url_only() -> None:
    text = dockerfile_text()

    if re.search(r"\bARG\s+\S+_SHA256\b", text):
        fail("Dockerfile downloads must not require user-maintained SHA256 build args")
    if "sha256sum -c -" in text or "download_verify" in text:
        fail("Dockerfile downloads must not be checksum-gated")

    expected_url_args = {
        "NCBI_BLAST_URL",
        "EDIRECT_TAR_URL",
        "EDIRECT_XTRACT_LINUX_URL",
        "EDIRECT_RCHIVE_LINUX_URL",
        "EDIRECT_TRANSMUTE_LINUX_URL",
        "WAIT_FOR_URL",
        "TRIMAL_URL",
        "MVIEW_URL",
        "RPSBPROC_URL",
        "CDDID_TBL_URL",
        "CDTRACK_URL",
        "FAMILY_SUPERFAMILY_LINKS_URL",
        "CDDANNOT_URL",
        "CDDANNOT_GENERIC_URL",
        "BITSCORE_SPECIFIC_URL",
        "MINICONDA_URL",
        "MICROMAMBA_URL",
    }
    missing = sorted(name for name in expected_url_args if name not in text)
    if missing:
        fail(f"Dockerfile is missing expected URL build args: {missing}")


def assert_non_root_runtime() -> None:
    app = read(ROOT / "Dockerfile")
    nginx = read(ROOT / "nginx" / "Dockerfile")
    production = read(ROOT / "docker-compose-production.yml")
    nginx_conf = read(ROOT / "nginx" / "nginx.conf")

    if "USER 10001:10001" not in app:
        fail("Application image must run as non-root UID/GID 10001")
    if "USER nginx" not in nginx:
        fail("Nginx image must run as the nginx user")
    if "--allow-root" in production:
        fail("Jupyter command must not use --allow-root")
    if "listen 8080;" not in nginx_conf:
        fail("Nginx must listen on an unprivileged container port")
    if "1337:8080" not in production:
        fail("Production Compose must publish host port 1337 to Nginx container port 8080")


def main() -> int:
    try:
        assert_base_images_are_explicit()
        assert_no_remote_shell_execution()
        assert_remote_downloads_are_url_only()
        assert_non_root_runtime()
    except AssertionError as exc:
        print(f"Docker build hardening check failed: {exc}", file=sys.stderr)
        return 1
    print("Docker build hardening check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
