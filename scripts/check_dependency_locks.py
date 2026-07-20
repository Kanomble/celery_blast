"""Validate CATHI dependency lock policy without contacting package indexes."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def fail(message: str) -> None:
    raise AssertionError(message)


def requirement_names(path: str) -> dict[str, str]:
    requirements: dict[str, str] = {}
    for line in read(path).splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        match = re.match(r"^([A-Za-z0-9_.-]+)==([^#\s]+)$", stripped)
        if not match:
            fail(f"{path} contains a non-exact requirement: {stripped}")
        requirements[match.group(1).lower().replace("_", "-")] = match.group(2)
    return requirements


def conda_versions(path: str) -> dict[str, str]:
    versions: dict[str, str] = {}
    in_dependencies = False
    for line in read(path).splitlines():
        stripped = line.strip()
        if stripped == "dependencies:":
            in_dependencies = True
            continue
        if not in_dependencies:
            continue
        if stripped and not line.startswith(" "):
            in_dependencies = False
            continue
        if not stripped.startswith("- "):
            continue
        spec = stripped[2:].strip()
        match = re.match(r"^([A-Za-z0-9_.-]+)=([^=\s]+)(?:=.*)?$", spec)
        if not match:
            fail(f"{path} contains a non-exact Conda dependency: {spec}")
        versions[match.group(1).lower().replace("_", "-")] = match.group(2)
    return versions


def assert_pip_lock() -> None:
    direct = requirement_names("requirements.in")
    locked = requirement_names("requirements.lock.txt")

    for name, version in direct.items():
        if locked.get(name) != version:
            fail(f"requirements.lock.txt does not lock direct dependency {name}=={version}")

    if "wget" in direct or "wget" in locked:
        fail("The unused Python wget package must not be in pip requirements")

    compatibility_entrypoint = read("requirements.txt").strip().splitlines()
    if compatibility_entrypoint != [
        "# Compatibility entry point for existing install commands.",
        "# Edit requirements.in, then regenerate and review requirements.lock.txt.",
        "-r requirements.lock.txt",
    ]:
        fail("requirements.txt must delegate to requirements.lock.txt")


def assert_conda_lock() -> None:
    for path in ("environment.yml", "environment-linux-64.lock.yml"):
        if "channel_priority:" in read(path):
            fail(f"{path} must not contain channel_priority; configure it in Conda")

    reviewed = conda_versions("environment.yml")
    locked = conda_versions("environment-linux-64.lock.yml")
    required = {
        "python": "3.8.18",
        "biopython": "1.78",
        "pandas": "2.0.3",
        "matplotlib": "3.7.5",
        "seaborn": "0.12.0",
        "bokeh": "2.4.3",
        "pyyaml": "6.0.2",
        "snakemake": "7.32.4",
        "notebook": "7.1.0",
        "pyopenssl": "24.0.0",
        "cryptography": "42.0.5",
        "aioeasywebdav": "2.4.0",
        "setuptools-scm": "8.0.4",
    }

    for name, version in required.items():
        if reviewed.get(name) != version:
            fail(f"environment.yml must pin {name}={version}")
        if locked.get(name) != version:
            fail(f"environment-linux-64.lock.yml must pin {name}={version}")

    if locked["cryptography"].startswith("38."):
        fail("cryptography 38.x is incompatible with pyOpenSSL 24.0.0")


def assert_dockerfile_policy() -> None:
    dockerfile = read("Dockerfile")
    forbidden = [
        "conda update --all",
        "pip install cryptography==38.0.4",
        "python3 -m pip install cryptography==38.0.4",
    ]
    for text in forbidden:
        if text in dockerfile:
            fail(f"Dockerfile still contains forbidden dependency command: {text}")

    ad_hoc_conda_installs = [
        "conda install -c bioconda snakemake",
        "conda install -c conda-forge notebook",
    ]
    for text in ad_hoc_conda_installs:
        if text in dockerfile:
            fail(f"Dockerfile must install runtime Conda packages from environment.yml, not {text}")
    if "python3 -m pip install --no-cache-dir -r requirements.lock.txt" not in dockerfile:
        fail("Dockerfile must install reviewed Python lock file")
    if "python3 -m pip check" not in dockerfile:
        fail("Dockerfile must run pip check after installing Python dependencies")
    if "micromamba create -y -p \"${CATHI_CONDA_ENV}\"" not in dockerfile:
        fail("Dockerfile must install Conda workflow tools into an isolated prefix")
    workflow_imports = ("from Bio import Entrez", "pandas", "matplotlib", "seaborn", "bokeh", "yaml")
    for import_text in workflow_imports:
        if import_text not in dockerfile:
            fail(f"Dockerfile must verify {import_text} is importable from the Conda workflow prefix")
    if "/blast/conda-envs/cathi-runtime" not in dockerfile:
        fail("Dockerfile must document the Conda workflow tool prefix")
    if "build_metadata/dependency_versions.txt" not in dockerfile:
        fail("Dockerfile must record dependency versions in build metadata")


def main() -> int:
    try:
        assert_pip_lock()
        assert_conda_lock()
        assert_dockerfile_policy()
    except AssertionError as exc:
        print(f"dependency lock check failed: {exc}", file=sys.stderr)
        return 1
    print("dependency lock check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
