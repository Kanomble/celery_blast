from blast_project.py_services import check_if_taxdb_exists
from django.test import TestCase

class PyServicesTestCase(TestCase):
    def test_check_if_taxdb_exists(self):
        bool = check_if_taxdb_exists()
        self.assertEqual(True,bool)