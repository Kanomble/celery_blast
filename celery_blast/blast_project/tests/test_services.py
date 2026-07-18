from pathlib import Path
from tempfile import TemporaryDirectory

from blast_project.py_services import cdd_database_prefix_exists, check_if_taxdb_exists
from django.test import SimpleTestCase, TestCase, tag


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


@tag('biological')
class PyServicesTestCase(TestCase):
    def test_check_if_taxdb_exists(self):
        bool = check_if_taxdb_exists()
        self.assertEqual(True,bool)
