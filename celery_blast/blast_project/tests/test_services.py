from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from blast_project.py_services import cdd_database_prefix_exists, check_if_taxdb_exists, delete_domain_database
from django.test import SimpleTestCase, TestCase, tag
from django_celery_results.models import TaskResult
from external_tools.models import DomainDatabase


class CddDatabaseStatusHelperTests(SimpleTestCase):
    def test_cdd_database_prefix_exists_detects_flat_blast_database_files(self):
        with TemporaryDirectory() as tmpdir:
            cdd_root = Path(tmpdir) / "CDD"
            cdd_root.mkdir()
            (cdd_root / "Cdd.pal").write_text("TITLE Cdd\n", encoding="utf-8")

            self.assertTrue(cdd_database_prefix_exists(str(cdd_root)))

    def test_cdd_database_prefix_exists_rejects_missing_prefix_files(self):
        with TemporaryDirectory() as tmpdir:
            cdd_root = Path(tmpdir) / "CDD"
            cdd_root.mkdir()
            (cdd_root / "README").write_text("not a BLAST database", encoding="utf-8")

            self.assertFalse(cdd_database_prefix_exists(str(cdd_root)))


class DomainDatabaseDeletionTests(TestCase):
    def test_delete_domain_database_removes_nested_cdd_directory_and_resets_state(self):
        with TemporaryDirectory() as tmpdir:
            database_root = Path(tmpdir) / "databases"
            cdd_root = database_root / "CDD"
            nested_root = cdd_root / "nested"
            nested_root.mkdir(parents=True)
            (nested_root / "Cdd.pal").write_text("TITLE Cdd\n", encoding="utf-8")
            task_result = TaskResult.objects.create(task_id="failed-cdd-task", status="FAILURE")
            DomainDatabase.objects.create(
                domain_database_loaded=False,
                domain_database_download_task_result=task_result,
            )

            with patch("blast_project.py_services.BLAST_DATABASE_DIR", str(database_root) + "/"):
                self.assertEqual(0, delete_domain_database())

            self.assertFalse(cdd_root.exists())
            domain_database = DomainDatabase.objects.get()
            self.assertFalse(domain_database.domain_database_loaded)
            self.assertIsNone(domain_database.domain_database_download_task_result)

    def test_delete_domain_database_removes_stray_cdd_file_and_resets_state(self):
        with TemporaryDirectory() as tmpdir:
            database_root = Path(tmpdir) / "databases"
            database_root.mkdir()
            cdd_file = database_root / "CDD"
            cdd_file.write_text("partial failed download", encoding="utf-8")
            task_result = TaskResult.objects.create(task_id="failed-cdd-task", status="FAILURE")
            DomainDatabase.objects.create(
                domain_database_loaded=False,
                domain_database_download_task_result=task_result,
            )

            with patch("blast_project.py_services.BLAST_DATABASE_DIR", str(database_root) + "/"):
                self.assertEqual(0, delete_domain_database())

            self.assertFalse(cdd_file.exists())
            domain_database = DomainDatabase.objects.get()
            self.assertFalse(domain_database.domain_database_loaded)
            self.assertIsNone(domain_database.domain_database_download_task_result)


@tag('biological')
class PyServicesTestCase(TestCase):
    def test_check_if_taxdb_exists(self):
        bool = check_if_taxdb_exists()
        self.assertEqual(True,bool)
