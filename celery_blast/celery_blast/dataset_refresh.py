"""Atomic refresh helpers for shared scientific reference datasets."""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import tarfile
import tempfile
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional
from urllib.request import urlopen


class DatasetRefreshError(Exception):
    """Base class for dataset refresh failures."""


class DatasetRefreshLocked(DatasetRefreshError):
    """Raised when another refresh for the same dataset is active."""


class DatasetVerificationError(DatasetRefreshError):
    """Raised when a staged dataset fails validation."""


DownloadCallable = Callable[[str, Path], None]


@dataclass(frozen=True)
class DatasetRefreshSpec:
    name: str
    source_url: str
    public_root: Path
    required_files: tuple[str, ...]
    archive_name: str
    archive_type: str = "tar.gz"
    expose_as: str = "files"
    minimum_total_bytes: int = 1
    timeout: int = 800

    def __post_init__(self):
        if self.archive_type not in {"tar.gz", "file"}:
            raise ValueError("archive_type must be 'tar.gz' or 'file'")
        if self.expose_as not in {"files", "directory"}:
            raise ValueError("expose_as must be 'files' or 'directory'")
        if not self.required_files:
            raise ValueError("required_files must not be empty")


def download_url(source_url: str, output_path: Path, timeout: int = 800) -> None:
    with urlopen(source_url, timeout=timeout) as response, output_path.open("wb") as output:
        shutil.copyfileobj(response, output)


def _metadata_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _store_parent(spec: DatasetRefreshSpec) -> Path:
    if spec.expose_as == "directory":
        return spec.public_root.parent
    return spec.public_root


def _dataset_store(spec: DatasetRefreshSpec) -> Path:
    return _store_parent(spec) / ".cathi_datasets" / spec.name


def active_pointer(spec: DatasetRefreshSpec) -> Path:
    return _dataset_store(spec) / "active"


def active_metadata_path(spec: DatasetRefreshSpec) -> Path:
    return active_pointer(spec) / ".cathi-dataset.json"


@contextlib.contextmanager
def dataset_refresh_lock(spec: DatasetRefreshSpec):
    lock_dir = _store_parent(spec) / ".cathi_locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    lock_path = lock_dir / f"{spec.name}.lock"
    try:
        fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise DatasetRefreshLocked(f"{spec.name} refresh is already in progress") from exc

    try:
        with os.fdopen(fd, "w", encoding="utf-8") as lock_file:
            lock_file.write(json.dumps({"pid": os.getpid(), "created_at": _metadata_now()}))
        yield
    finally:
        with contextlib.suppress(FileNotFoundError):
            lock_path.unlink()


def _safe_tar_members(archive_path: Path) -> Iterable[tarfile.TarInfo]:
    try:
        with tarfile.open(archive_path, "r:*") as archive:
            for member in archive.getmembers():
                member_path = Path(member.name)
                if member_path.is_absolute() or ".." in member_path.parts:
                    raise DatasetVerificationError(f"archive contains unsafe path: {member.name}")
                yield member
    except tarfile.TarError as exc:
        raise DatasetVerificationError(f"archive is not readable: {archive_path.name}") from exc


def _extract_tar_safely(archive_path: Path, destination: Path) -> None:
    members = list(_safe_tar_members(archive_path))
    with tarfile.open(archive_path, "r:*") as archive:
        for member in members:
            target = (destination / member.name).resolve()
            if destination.resolve() not in (target, *target.parents):
                raise DatasetVerificationError(f"archive member escapes destination: {member.name}")
        archive.extractall(destination, members=members)


def _verify_required_files(extracted_root: Path, required_files: Iterable[str], minimum_total_bytes: int) -> None:
    total_size = 0
    for required in required_files:
        path = extracted_root / required
        if not path.exists():
            raise DatasetVerificationError(f"required dataset path is missing: {required}")
        if path.is_file():
            size = path.stat().st_size
            if size == 0:
                raise DatasetVerificationError(f"required dataset file is empty: {required}")
            total_size += size
        elif path.is_dir():
            files = [item for item in path.rglob("*") if item.is_file()]
            if not files:
                raise DatasetVerificationError(f"required dataset directory is empty: {required}")
            total_size += sum(item.stat().st_size for item in files)
    if total_size < minimum_total_bytes:
        raise DatasetVerificationError("dataset did not meet minimum size sanity check")


def _replace_symlink(link_path: Path, target: Path) -> None:
    tmp_link = link_path.with_name(f".{link_path.name}.tmp-{uuid.uuid4().hex}")
    if target.is_absolute():
        link_target = target
    else:
        link_target = os.path.relpath(target, start=link_path.parent)
    os.symlink(str(link_target), str(tmp_link), target_is_directory=target.is_dir())
    os.replace(str(tmp_link), str(link_path))


def _copy_legacy_public_files(spec: DatasetRefreshSpec, legacy_root: Path, version_dir: Path) -> None:
    if spec.expose_as == "directory":
        if spec.public_root.exists() and not spec.public_root.is_symlink():
            legacy_target = version_dir.parent / f"legacy-{uuid.uuid4().hex}"
            shutil.copytree(spec.public_root, legacy_target)
        return

    for required in spec.required_files:
        public_file = legacy_root / required
        if public_file.exists() and not public_file.is_symlink():
            legacy_file = version_dir.parent / f"legacy-{uuid.uuid4().hex}" / required
            legacy_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(public_file, legacy_file)


def _expose_active_version(spec: DatasetRefreshSpec) -> None:
    active = active_pointer(spec)
    if spec.expose_as == "directory":
        if spec.public_root.exists() and not spec.public_root.is_symlink():
            legacy_target = _dataset_store(spec) / "versions" / f"legacy-public-{uuid.uuid4().hex}"
            os.replace(str(spec.public_root), str(legacy_target))
        _replace_symlink(spec.public_root, active)
        return

    spec.public_root.mkdir(parents=True, exist_ok=True)
    for required in spec.required_files:
        public_file = spec.public_root / required
        public_file.parent.mkdir(parents=True, exist_ok=True)
        _replace_symlink(public_file, active / required)

    metadata_link = spec.public_root / f"{spec.name}.dataset.json"
    _replace_symlink(metadata_link, active / ".cathi-dataset.json")


def refresh_dataset(spec: DatasetRefreshSpec, download: Optional[DownloadCallable] = None) -> dict:
    parent = _store_parent(spec)
    parent.mkdir(parents=True, exist_ok=True)
    store = _dataset_store(spec)
    versions_dir = store / "versions"
    versions_dir.mkdir(parents=True, exist_ok=True)
    staging_parent = parent / ".cathi_refresh"
    staging_parent.mkdir(parents=True, exist_ok=True)

    download_func = download or (lambda url, path: download_url(url, path, timeout=spec.timeout))
    staging_dir = Path(tempfile.mkdtemp(prefix=f"{spec.name}-", dir=str(staging_parent)))
    download_started = _metadata_now()

    try:
        with dataset_refresh_lock(spec):
            archive_path = staging_dir / spec.archive_name
            extracted_root = staging_dir / "extracted"
            extracted_root.mkdir()

            download_func(spec.source_url, archive_path)

            if spec.archive_type == "tar.gz":
                _extract_tar_safely(archive_path, extracted_root)
            else:
                target = extracted_root / spec.required_files[0]
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(archive_path, target)

            _verify_required_files(extracted_root, spec.required_files, spec.minimum_total_bytes)

            version_id = f"{int(time.time())}-{uuid.uuid4().hex}"
            version_dir = versions_dir / version_id
            metadata = {
                "dataset": spec.name,
                "version": version_id,
                "source": spec.source_url,
                "download_started_at": download_started,
                "download_completed_at": _metadata_now(),
                "activated_at": None,
                "required_files": list(spec.required_files),
            }
            (extracted_root / ".cathi-dataset.json").write_text(
                json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            _copy_legacy_public_files(spec, parent, version_dir)
            os.replace(str(extracted_root), str(version_dir))
            metadata["activated_at"] = _metadata_now()
            (version_dir / ".cathi-dataset.json").write_text(
                json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            active_tmp = store / f"active.tmp-{uuid.uuid4().hex}"
            os.symlink(os.path.relpath(version_dir, start=store), active_tmp, target_is_directory=True)
            os.replace(str(active_tmp), str(active_pointer(spec)))
            _expose_active_version(spec)
            return metadata
    finally:
        shutil.rmtree(staging_dir, ignore_errors=True)
