"""Static checks for Docker build hardening policy."""

from __future__ import annotations

import json
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


def load_manifest() -> dict[str, object]:
    return json.loads((ROOT / "build" / "remote-artifacts.json").read_text(encoding="utf-8"))


def assert_base_images_are_digest_pinned() -> None:
    for path in DOCKERFILES:
        for line in read(path).splitlines():
            stripped = line.strip()
            if not stripped.startswith("FROM "):
                continue
            if "${UBUNTU_FOCAL_AMD64_IMAGE}" in stripped:
                continue
            if "@sha256:" not in stripped:
                fail(f"{path} has an unpinned base image: {stripped}")

    dockerfile = read(ROOT / "Dockerfile")
    if "UBUNTU_FOCAL_AMD64_IMAGE=ubuntu:focal@sha256:" not in dockerfile:
        fail("Dockerfile must pin ubuntu:focal by digest through UBUNTU_FOCAL_AMD64_IMAGE")


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


def assert_artifacts_are_verified() -> None:
    manifest = load_manifest()
    text = dockerfile_text()

    for artifact in manifest["artifacts"]:
        artifact_id = artifact["id"]
        url = artifact["url"]
        sha256 = artifact["sha256"]
        if not str(url).startswith("https://"):
            fail(f"{artifact_id} must use HTTPS: {url}")
        if not re.fullmatch(r"[0-9a-f]{64}", sha256):
            fail(f"{artifact_id} must have a SHA256 checksum")
        if artifact_id == "miniconda":
            expected_names = ["MINICONDA_URL", "MINICONDA_SHA256"]
        elif artifact_id == "edirect":
            expected_names = ["EDIRECT_TAR_URL", "EDIRECT_TAR_SHA256"]
        else:
            expected_names = [
                artifact_id.upper().replace("-", "_").replace(".", "_") + "_URL",
                artifact_id.upper().replace("-", "_").replace(".", "_") + "_SHA256",
            ]
        for value in (url, sha256):
            if value not in text:
                fail(f"{artifact_id} is missing from Dockerfile: {value}")
        if not any(name in text for name in expected_names):
            fail(f"{artifact_id} must be represented by checksum build args")

    if "sha256sum -c -" not in text:
        fail("Dockerfile must fail builds through sha256sum -c -")


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


def assert_sbom_includes_principal_tools() -> None:
    manifest = load_manifest()
    principal = {
        artifact["name"]
        for artifact in manifest["artifacts"]
        if artifact.get("principal_scientific_tool")
    }
    expected = {"NCBI BLAST+", "NCBI EDirect", "trimAl", "MView", "RpsbProc"}
    missing = expected - principal
    if missing:
        fail(f"remote artifact manifest is missing principal tools: {sorted(missing)}")

    sbom_script = read(ROOT / "scripts" / "generate_release_sbom.py")
    dockerfile = read(ROOT / "Dockerfile")
    if "CycloneDX" not in sbom_script:
        fail("SBOM generator must produce CycloneDX output")
    if "scientific-tools.sbom.cdx.json" not in dockerfile:
        fail("Dockerfile must generate the scientific tools SBOM")


def main() -> int:
    try:
        assert_base_images_are_digest_pinned()
        assert_no_remote_shell_execution()
        assert_artifacts_are_verified()
        assert_non_root_runtime()
        assert_sbom_includes_principal_tools()
    except AssertionError as exc:
        print(f"Docker build hardening check failed: {exc}", file=sys.stderr)
        return 1
    print("Docker build hardening check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
