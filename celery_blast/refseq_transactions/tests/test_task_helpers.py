from tempfile import TemporaryDirectory
from unittest.mock import ANY, call, patch

from django.test import SimpleTestCase
from django.test.utils import override_settings

from refseq_transactions.tasks import (
    build_download_assembly_summary_command,
    build_download_and_gunzip_command,
    build_format_blast_database_command,
    download_wget_ftp_paths,
    format_blast_databases,
    run_download_assembly_summary_command,
    run_download_and_gunzip_command,
    run_format_blast_database_command,
)


class FakeProgressRecorder:
    def __init__(self):
        self.progress_updates = []

    def set_progress(self, current, total, description):
        self.progress_updates.append((current, total, description))


class RefseqTaskHelperTests(SimpleTestCase):
    def test_build_download_and_gunzip_command_uses_argument_passing(self):
        command = build_download_and_gunzip_command(
            'https://example.test/protein.faa.gz',
            '/data/protein.faa',
        )

        self.assertEqual(
            [
                'bash',
                '-o', 'pipefail',
                '-c', 'wget -qO- "$1" | gzip -d > "$2"',
                'download_assembly',
                'https://example.test/protein.faa.gz',
                '/data/protein.faa',
            ],
            command,
        )

    def test_run_download_and_gunzip_command_delegates_to_runner_with_checking(self):
        class Result:
            returncode = 0

        command = build_download_and_gunzip_command('url', '/data/out.faa')
        with patch('refseq_transactions.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_download_and_gunzip_command(command)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(
            command,
            timeout=300,
            shell=False,
            logger=ANY,
            check=True,
            cleanup_exceptions=ANY,
        )

    def test_build_download_assembly_summary_command(self):
        command = build_download_assembly_summary_command(
            'https://example.test/assembly_summary_refseq.txt',
            '/data/assembly_summary_refseq.txt',
        )

        self.assertEqual(
            [
                'curl',
                '--fail',
                '--location',
                '--show-error',
                '--silent',
                '-o',
                '/data/assembly_summary_refseq.txt',
                'https://example.test/assembly_summary_refseq.txt',
            ],
            command,
        )

    def test_run_download_assembly_summary_command_delegates_to_runner_with_checking(self):
        class Result:
            returncode = 0

        command = build_download_assembly_summary_command('url', '/data/summary.txt')
        with patch('refseq_transactions.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_download_assembly_summary_command(command, timeout=45)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(
            command,
            timeout=45,
            shell=False,
            logger=ANY,
            check=True,
            cleanup_exceptions=ANY,
        )

    def test_build_format_blast_database_command(self):
        command = build_format_blast_database_command(
            '/data/database_chunk_1.faa',
            '/data/acc_taxmap_file_1.table',
        )

        self.assertEqual(
            [
                'makeblastdb',
                '-in', '/data/database_chunk_1.faa',
                '-dbtype', 'prot',
                '-taxid_map', '/data/acc_taxmap_file_1.table',
                '-parse_seqids',
                '-out', '/data/database_chunk_1.faa',
                '-blastdb_version', '5',
            ],
            command,
        )

    @override_settings(SUBPROCESS_TIME_LIMIT=75)
    def test_run_format_blast_database_command_delegates_to_runner_without_checking_returncode(self):
        class Result:
            returncode = 2

        command = ['makeblastdb', '-in', 'db.faa']
        with patch('refseq_transactions.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_format_blast_database_command(command)

        self.assertEqual(2, returncode)
        run_command.assert_called_once_with(
            command,
            timeout=75,
            shell=False,
            logger=ANY,
            check=False,
            cleanup_exceptions=ANY,
        )

    def test_format_blast_databases_tracks_successful_and_failed_chunks(self):
        progress_recorder = FakeProgressRecorder()

        with patch('refseq_transactions.tasks.run_format_blast_database_command',
                   side_effect=[0, 1]) as run_command:
            available_databases, errorlist = format_blast_databases(
                '/data/',
                [1, 2],
                progress_recorder,
            )

        self.assertEqual(['/data/database_chunk_1.faa'], available_databases)
        self.assertEqual([2], errorlist)
        self.assertEqual(
            [
                call([
                    'makeblastdb',
                    '-in', '/data/database_chunk_1.faa',
                    '-dbtype', 'prot',
                    '-taxid_map', '/data/acc_taxmap_file_1.table',
                    '-parse_seqids',
                    '-out', '/data/database_chunk_1.faa',
                    '-blastdb_version', '5',
                ]),
                call([
                    'makeblastdb',
                    '-in', '/data/database_chunk_2.faa',
                    '-dbtype', 'prot',
                    '-taxid_map', '/data/acc_taxmap_file_2.table',
                    '-parse_seqids',
                    '-out', '/data/database_chunk_2.faa',
                    '-blastdb_version', '5',
                ]),
            ],
            run_command.call_args_list,
        )
        self.assertEqual((75, 100, 'starting to format database chunks'), progress_recorder.progress_updates[0])
        self.assertEqual((99, 100, 'formatted chunk 2'), progress_recorder.progress_updates[-1])

    def test_download_wget_ftp_paths_uses_safe_download_command(self):
        progress_recorder = FakeProgressRecorder()
        ftp_path = 'https://example.test/GCF_000001_protein.faa.gz'

        with TemporaryDirectory() as tmpdir, \
                patch('refseq_transactions.tasks.run_download_and_gunzip_command', return_value=0) as run_command, \
                patch('refseq_transactions.tasks.chmod') as chmod:
            downloaded_files, errorlist = download_wget_ftp_paths(
                tmpdir + '/',
                {ftp_path: '1140'},
                progress_recorder,
            )

        self.assertEqual({'GCF_000001_protein.faa': '1140'}, downloaded_files)
        self.assertEqual([], errorlist)
        run_command.assert_called_once_with(
            [
                'bash',
                '-o', 'pipefail',
                '-c', 'wget -qO- "$1" | gzip -d > "$2"',
                'download_assembly',
                ftp_path,
                ANY,
            ]
        )
        self.assertTrue(run_command.call_args[0][0][-1].endswith('/GCF_000001_protein.faa'))
        chmod.assert_called_once()
