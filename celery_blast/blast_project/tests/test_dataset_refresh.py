import os
import shutil
import tarfile
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

from django.test import SimpleTestCase

from celery_blast.dataset_refresh import (
    DatasetRefreshLocked,
    DatasetRefreshSpec,
    DatasetVerificationError,
    active_pointer,
    dataset_refresh_lock,
    refresh_dataset,
    sha256_file,
)


def symlink_supported():
    with tempfile.TemporaryDirectory() as tmpdir:
        target = Path(tmpdir) / "target"
        link = Path(tmpdir) / "link"
        target.write_text("ok", encoding="utf-8")
        try:
            os.symlink(target, link)
        except (OSError, NotImplementedError):
            return False
        return link.read_text(encoding="utf-8") == "ok"


def make_tar(path, entries):
    with tarfile.open(path, "w:gz") as archive:
        for name, data in entries.items():
            source = Path(path).parent / name.replace("/", "_")
            source.write_bytes(data)
            archive.add(source, arcname=name)


def copy_download(source, destination):
    shutil.copy2(source, destination)


@unittest.skipUnless(symlink_supported(), "atomic dataset activation requires symlink support")
class DatasetRefreshTests(SimpleTestCase):
    def make_spec(self, root, archive, checksum=None):
        return DatasetRefreshSpec(
            name="taxdb",
            source_url=str(archive),
            public_root=Path(root) / "databases",
            required_files=("taxdb.btd", "taxdb.bti"),
            expected_sha256=checksum or sha256_file(Path(archive)),
            archive_name="taxdb.tar.gz",
            archive_type="tar.gz",
            expose_as="files",
            minimum_total_bytes=2,
        )

    def test_successful_refresh_activates_new_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            archive = Path(tmpdir) / "taxdb.tar.gz"
            make_tar(archive, {"taxdb.btd": b"new-btd", "taxdb.bti": b"new-bti"})
            spec = self.make_spec(tmpdir, archive)

            metadata = refresh_dataset(spec, download=copy_download)

            self.assertEqual("taxdb", metadata["dataset"])
            self.assertEqual("new-btd", (spec.public_root / "taxdb.btd").read_text(encoding="utf-8"))
            self.assertTrue((active_pointer(spec) / ".cathi-dataset.json").is_file())

    def test_interrupted_download_leaves_old_version_active(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            old_archive = Path(tmpdir) / "old.tar.gz"
            new_archive = Path(tmpdir) / "new.tar.gz"
            make_tar(old_archive, {"taxdb.btd": b"old-btd", "taxdb.bti": b"old-bti"})
            make_tar(new_archive, {"taxdb.btd": b"new-btd", "taxdb.bti": b"new-bti"})
            refresh_dataset(self.make_spec(tmpdir, old_archive), download=copy_download)

            def interrupted_download(source, destination):
                destination.write_bytes(b"partial")
                raise RuntimeError("simulated interruption")

            with self.assertRaises(RuntimeError):
                refresh_dataset(self.make_spec(tmpdir, new_archive), download=interrupted_download)

            self.assertEqual("old-btd", (Path(tmpdir) / "databases" / "taxdb.btd").read_text(encoding="utf-8"))

    def test_checksum_failure_leaves_old_version_active(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            old_archive = Path(tmpdir) / "old.tar.gz"
            new_archive = Path(tmpdir) / "new.tar.gz"
            make_tar(old_archive, {"taxdb.btd": b"old-btd", "taxdb.bti": b"old-bti"})
            make_tar(new_archive, {"taxdb.btd": b"new-btd", "taxdb.bti": b"new-bti"})
            refresh_dataset(self.make_spec(tmpdir, old_archive), download=copy_download)
            bad_spec = self.make_spec(tmpdir, new_archive, checksum="0" * 64)

            with self.assertRaises(DatasetVerificationError):
                refresh_dataset(bad_spec, download=copy_download)

            self.assertEqual("old-btd", (Path(tmpdir) / "databases" / "taxdb.btd").read_text(encoding="utf-8"))

    def test_invalid_archive_contents_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            old_archive = Path(tmpdir) / "old.tar.gz"
            bad_archive = Path(tmpdir) / "bad.tar.gz"
            make_tar(old_archive, {"taxdb.btd": b"old-btd", "taxdb.bti": b"old-bti"})
            make_tar(bad_archive, {"taxdb.btd": b"only-one-file"})
            refresh_dataset(self.make_spec(tmpdir, old_archive), download=copy_download)

            with self.assertRaises(DatasetVerificationError):
                refresh_dataset(self.make_spec(tmpdir, bad_archive), download=copy_download)

            self.assertEqual("old-btd", (Path(tmpdir) / "databases" / "taxdb.btd").read_text(encoding="utf-8"))

    def test_archive_traversal_entries_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            old_archive = Path(tmpdir) / "old.tar.gz"
            traversal_archive = Path(tmpdir) / "traversal.tar.gz"
            make_tar(old_archive, {"taxdb.btd": b"old-btd", "taxdb.bti": b"old-bti"})
            make_tar(traversal_archive, {"../escape": b"bad", "taxdb.btd": b"x", "taxdb.bti": b"y"})
            refresh_dataset(self.make_spec(tmpdir, old_archive), download=copy_download)

            with self.assertRaises(DatasetVerificationError):
                refresh_dataset(self.make_spec(tmpdir, traversal_archive), download=copy_download)

            self.assertFalse((Path(tmpdir) / "escape").exists())
            self.assertEqual("old-btd", (Path(tmpdir) / "databases" / "taxdb.btd").read_text(encoding="utf-8"))

    def test_concurrent_refresh_attempts_do_not_overlap(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            archive = Path(tmpdir) / "taxdb.tar.gz"
            make_tar(archive, {"taxdb.btd": b"new-btd", "taxdb.bti": b"new-bti"})
            spec = self.make_spec(tmpdir, archive)

            with dataset_refresh_lock(spec):
                with self.assertRaises(DatasetRefreshLocked):
                    refresh_dataset(spec, download=copy_download)

    def test_activation_uses_atomic_filesystem_operation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            archive = Path(tmpdir) / "taxdb.tar.gz"
            make_tar(archive, {"taxdb.btd": b"new-btd", "taxdb.bti": b"new-bti"})
            spec = self.make_spec(tmpdir, archive)

            with patch("celery_blast.dataset_refresh.os.replace", wraps=os.replace) as replace:
                refresh_dataset(spec, download=copy_download)

            active = str(active_pointer(spec))
            self.assertTrue(any(call.args[1] == active for call in replace.call_args_list))

    def test_directory_dataset_refresh_retains_legacy_public_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            public = Path(tmpdir) / "databases" / "CDD"
            public.mkdir(parents=True)
            (public / "old.txt").write_text("old-cdd", encoding="utf-8")
            archive = Path(tmpdir) / "cdd.tar.gz"
            make_tar(archive, {"Cdd/Cdd.pal": b"new-cdd"})
            spec = DatasetRefreshSpec(
                name="cdd",
                source_url=str(archive),
                public_root=public,
                required_files=("Cdd",),
                expected_sha256=sha256_file(archive),
                archive_name="Cdd_LE.tar.gz",
                archive_type="tar.gz",
                expose_as="directory",
                minimum_total_bytes=1,
            )

            refresh_dataset(spec, download=copy_download)

            self.assertEqual("new-cdd", (public / "Cdd" / "Cdd.pal").read_text(encoding="utf-8"))
            legacy_dirs = list((Path(tmpdir) / "databases" / ".cathi_datasets" / "cdd" / "versions").glob("legacy-*"))
            self.assertTrue(any((legacy / "old.txt").read_text(encoding="utf-8") == "old-cdd" for legacy in legacy_dirs))

    def test_unreadable_archive_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            archive = Path(tmpdir) / "not-a-tar.gz"
            archive.write_bytes(b"not a tar")
            spec = self.make_spec(tmpdir, archive)

            with self.assertRaises(DatasetVerificationError):
                refresh_dataset(spec, download=copy_download)
