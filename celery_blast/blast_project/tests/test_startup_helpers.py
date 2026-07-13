from unittest.mock import patch

from django.test import SimpleTestCase

from scripts.startup import (
    build_startup_taxdb_download_command,
    build_startup_taxdb_extract_command,
    run_startup_command,
)


class StartupHelperTests(SimpleTestCase):
    def test_build_startup_taxdb_download_command(self):
        command = build_startup_taxdb_download_command(
            'https://example.test/taxdb.tar.gz',
            '/data/taxdb.tar.gz',
        )

        self.assertEqual(
            ['wget', 'https://example.test/taxdb.tar.gz', '-q', '-O', '/data/taxdb.tar.gz'],
            command,
        )

    def test_build_startup_taxdb_extract_command(self):
        command = build_startup_taxdb_extract_command('/data/taxdb.tar.gz', '/data/databases/')

        self.assertEqual(['tar', '-zxvf', '/data/taxdb.tar.gz', '-C', '/data/databases/'], command)

    def test_run_startup_command_delegates_to_shared_runner(self):
        class Result:
            returncode = 0

        command = ['wget', 'url']
        with patch('scripts.startup.run_external_command', return_value=Result()) as run_command:
            returncode = run_startup_command(command, 600)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(command, timeout=600, shell=False, check=True)
