import importlib.util
import shutil
import subprocess
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path
from unittest.mock import ANY, mock_open, patch

from django.contrib.auth.models import AnonymousUser
from django.template.loader import render_to_string
from django.test import RequestFactory, SimpleTestCase

from external_tools.shiptv_html import make_tree_drawable, patch_shiptv_html_file
from external_tools.genbank_download_clinker_synteny import (
    build_download_genbank_command,
    run_genbank_download_command,
)
from external_tools.tasks import (
    build_clinker_command,
    build_cdd_download_command,
    build_cdd_extract_command,
    build_fasttree_command,
    build_mafft_command,
    build_rpsblast_command,
    build_shiptv_command,
    run_external_tool_command,
)
from external_tools.entrez_search_service import (
    build_efetch_fasta_from_input_command,
    build_entrez_search_command,
    build_protein_query_to_fasta_command,
    build_pubmed_to_fasta_command,
    run_edirect_command,
)


class SyntenyQuerySequenceStub:
    def __init__(self, accession_id):
        self.query_accession_id = accession_id
        self.query_accession_information = ''

    def check_if_synteny_calculation_task_is_complete(self):
        return 'PROGRESS'


class SyntenyQuerySequenceCollectionStub:
    def __init__(self, query_sequences):
        self._query_sequences = query_sequences

    def get_queryset(self):
        return self._query_sequences


class SyntenyExternalToolsStub:
    def __init__(self, query_sequences):
        self.query_sequences = SyntenyQuerySequenceCollectionStub(query_sequences)


class SyntenyDetectionDashboardTemplateTests(SimpleTestCase):
    def test_progress_controls_use_query_accession_id_consistently(self):
        query_sequences = [
            SyntenyQuerySequenceStub('WP_000001.1'),
            SyntenyQuerySequenceStub('WP_000002.1'),
        ]
        request = RequestFactory().get('/external_tools/synteny/')
        request.user = AnonymousUser()

        html = render_to_string(
            'external_tools/synteny_detection_dashboard.html',
            {
                'project_id': 7,
                'remote_or_local': 'local',
                'qseqids': SyntenyExternalToolsStub(query_sequences),
            },
            request=request,
        )

        self.assertNotIn('query_accession_Id', html)
        for query_sequence in query_sequences:
            accession_id = query_sequence.query_accession_id
            expected_ids = [
                f"synteny-calculation-progress-{accession_id}",
                f"synteny-calculation-button-{accession_id}",
                f"synteny-view-button-disabled-{accession_id}",
                f"synteny-view-button-{accession_id}",
                f"synteny-delete_button-disabled-{accession_id}",
                f"synteny-delete_button-{accession_id}",
            ]

            for expected_id in expected_ids:
                self.assertIn(expected_id, html)

            self.assertIn(
                f'document.getElementById("synteny-delete_button-disabled-{accession_id}")',
                html,
            )
            self.assertIn(
                f'document.getElementById("synteny-delete_button-{accession_id}")',
                html,
            )
            self.assertEqual(
                1,
                html.count(f'id="synteny-delete_button-disabled-{accession_id}"'),
            )
            self.assertEqual(
                1,
                html.count(f'id="synteny-delete_button-{accession_id}"'),
            )


class EntrezSearchServiceCommandTests(SimpleTestCase):
    def test_build_entrez_search_command_passes_query_as_argument(self):
        command = build_entrez_search_command(
            'protein',
            50,
            'kinase "quoted"',
            '/tmp/search.tsv',
        )

        self.assertEqual('bash', command[0])
        self.assertEqual('-o', command[1])
        self.assertEqual('pipefail', command[2])
        self.assertIn('esearch -db "$1" -query "$2"', command[4])
        self.assertEqual('protein', command[6])
        self.assertEqual('kinase "quoted"', command[7])
        self.assertEqual('50', command[8])
        self.assertEqual('Id Caption Title Organism', command[9])
        self.assertEqual('/tmp/search.tsv', command[10])

    def test_build_efetch_fasta_from_input_command_passes_paths_as_arguments(self):
        command = build_efetch_fasta_from_input_command('/tmp/ids.txt', '/tmp/out.faa')

        self.assertEqual(
            [
                'bash',
                '-c',
                'efetch -db protein -format fasta -input "$1" > "$2"',
                'efetch_fasta',
                '/tmp/ids.txt',
                '/tmp/out.faa',
            ],
            command,
        )

    def test_build_pubmed_to_fasta_command_preserves_active_download_shape(self):
        command = build_pubmed_to_fasta_command('hydra symbiont', '/tmp/out.faa')

        self.assertIn('elink -target protein | efetch -format fasta -start 1 -stop 100', command[4])
        self.assertEqual('hydra symbiont', command[6])
        self.assertEqual('/tmp/out.faa', command[7])

    def test_build_pubmed_to_fasta_command_can_use_refseq_filter(self):
        command = build_pubmed_to_fasta_command('hydra symbiont', '/tmp/out.faa', refseq_only=True)

        self.assertIn('elink -target protein | efilter -source refseq | efetch -format fasta', command[4])
        self.assertNotIn('-start 1 -stop 100', command[4])

    def test_build_protein_query_to_fasta_command(self):
        command = build_protein_query_to_fasta_command('WP_12345', '/tmp/out.faa')

        self.assertEqual('WP_12345', command[6])
        self.assertEqual('/tmp/out.faa', command[7])
        self.assertIn('esearch -db protein -query "$1"', command[4])

    def test_run_edirect_command_delegates_to_runner_without_checking_returncode(self):
        class Result:
            returncode = 2

        command = ['bash', '-c', 'exit 2']
        with patch('external_tools.entrez_search_service.run_external_command', return_value=Result()) as run_command:
            returncode = run_edirect_command(command, timeout=30)

        self.assertEqual(2, returncode)
        run_command.assert_called_once_with(command, timeout=30, shell=False, check=False)


class ExternalToolCommandTests(SimpleTestCase):
    def test_build_mafft_command_passes_paths_as_arguments(self):
        command = build_mafft_command('/tmp/input with spaces.faa')

        self.assertEqual(['mafft', '/tmp/input with spaces.faa'], command)

    def test_build_fasttree_command_passes_paths_as_arguments(self):
        command = build_fasttree_command('/tmp/input.msa')

        self.assertEqual(['fasttree', '-lg', '/tmp/input.msa'], command)

    def test_build_shiptv_command_uses_argv_without_shell(self):
        command = build_shiptv_command('/tmp/tree.nwk', '/tmp/metadata.tsv', '/tmp/tree.html')

        self.assertEqual(
            [
                'shiptv',
                '--newick', '/tmp/tree.nwk',
                '--metadata', '/tmp/metadata.tsv',
                '--output-html', '/tmp/tree.html',
            ],
            command,
        )

    def test_build_clinker_command_expands_to_argv(self):
        command = build_clinker_command(['/tmp/a.gbk', '/tmp/b.gbk'], '/tmp/clinker.html')

        self.assertEqual(
            [
                'clinker',
                '/tmp/a.gbk',
                '/tmp/b.gbk',
                '-p',
                '/tmp/clinker.html',
                '-i',
                '0.25',
            ],
            command,
        )

    def test_build_cdd_download_command(self):
        command = build_cdd_download_command('https://example.test/Cdd_LE.tar.gz', '/tmp/Cdd_LE.tar.gz')

        self.assertEqual(
            ['wget', 'https://example.test/Cdd_LE.tar.gz', '-q', '-O', '/tmp/Cdd_LE.tar.gz'],
            command,
        )

    def test_build_cdd_extract_command(self):
        command = build_cdd_extract_command('/tmp/Cdd_LE.tar.gz', '/tmp/CDD')

        self.assertEqual(['tar', '-zxvf', '/tmp/Cdd_LE.tar.gz', '-C', '/tmp/CDD/'], command)

    def test_build_rpsblast_command(self):
        command = build_rpsblast_command(
            '/tmp/query.faa',
            '/tmp/cdd',
            '/tmp/domains.tsv',
            {
                'rps_e_value': 0.001,
                'rps_num_threads': 4,
                'rps_max_hsps': 2,
                'rps_num_alignments': 100,
            },
        )

        self.assertEqual('rpsblast', command[0])
        self.assertEqual('/tmp/query.faa', command[2])
        self.assertEqual('/tmp/cdd', command[4])
        self.assertEqual('/tmp/domains.tsv', command[8])
        self.assertEqual('0.001', command[10])
        self.assertEqual('4', command[12])
        self.assertEqual('2', command[14])
        self.assertEqual('100', command[16])

    def test_run_external_tool_command_delegates_to_runner_with_check_enabled(self):
        class Result:
            returncode = 0

        command = ['shiptv', '--version']
        with patch('external_tools.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_external_tool_command(command, timeout=45)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(
            command,
            timeout=45,
            shell=False,
            logger=ANY,
            check=True,
            cleanup_exceptions=(),
        )

    def test_run_external_tool_command_redirects_stdout_to_path(self):
        class Result:
            returncode = 0

        command = ['mafft', '/tmp/input.faa']
        opener = mock_open()
        with patch('builtins.open', opener), \
                patch('external_tools.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_external_tool_command(command, stdout_path='/tmp/output.msa')

        self.assertEqual(0, returncode)
        opener.assert_called_once_with('/tmp/output.msa', 'w')
        run_command.assert_called_once_with(
            command,
            timeout=40000,
            shell=False,
            logger=ANY,
            check=True,
            cleanup_exceptions=(),
            stdout=opener.return_value.__enter__.return_value,
        )

    def test_run_external_tool_command_delegates_cleanup_exceptions(self):
        class Result:
            returncode = 0

        cleanup_exceptions = (RuntimeError,)
        command = ['rpsblast', '-query', '/tmp/query.faa']
        with patch('external_tools.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_external_tool_command(command, timeout=45, cleanup_exceptions=cleanup_exceptions)

        self.assertEqual(0, returncode)
        self.assertEqual(cleanup_exceptions, run_command.call_args.kwargs['cleanup_exceptions'])


class ShiptvHtmlPatchTests(SimpleTestCase):
    @unittest.skipUnless(shutil.which('shiptv'), 'shiptv CLI is not available')
    def test_generated_standalone_html_uses_public_ag_grid_api(self):
        html_path = self.generate_shiptv_html("(A:0.00000,B:0.00000):0.00000;")
        patch_shiptv_html_file(html_path)

        html = html_path.read_text(encoding='utf-8')
        app_code = self.extract_shiptv_app_code(html)

        self.assertIn('agGrid.createGrid', app_code)
        self.assertNotIn('new agGrid.Grid', app_code)
        self.assertNotIn('.context.beans', app_code)
        self.assertNotIn('.beanInstance', app_code)
        self.assertIn("mode: 'multiRow'", app_code)
        self.assertIn('enableClickSelection: true', app_code)
        self.assertIn('theme: agGrid.themeBalham', app_code)
        self.assertIn('var metadataColumns = [];', app_code)
        self.assertIn('metadataColumns.reduce', app_code)
        self.assertIn('node.setSelected(', app_code)
        self.assertIn("gridApi.setGridOption('rowData', newRowData);", app_code)
        self.assertIn('tree.load(makeTreeDrawable(newick_string));', app_code)

        self.assertIn('<meta charset="utf-8">', html)
        self.assertIn('<title>shiptv: Standalone HTML Interactive Phylogenetic Tree Visualization</title>', html)
        self.assertEqual(1, html.count('<head>'))
        self.assertIn('<div id="metadata-grid" style="height: 600px; width: 100%;"></div>', html)
        self.assertNotIn('class="ag-theme-balham" id="metadata-grid"', html)
        self.assertNotIn('overflow: none', html)
        self.assertEqual([], ExternalResourceParser.find_resources(html))
        self.assertNotIn('https://cdnjs.cloudflare.com/ajax/libs/chroma-js', html)

    def test_make_tree_drawable_converts_only_all_zero_branch_length_trees(self):
        self.assertEqual(
            '(A:1,B:1):1;',
            make_tree_drawable('(A:0.00000,B:0.00000):0.00000;'),
        )
        self.assertEqual(
            '(A:0.00000,B:0.5):0.00000;',
            make_tree_drawable('(A:0.00000,B:0.5):0.00000;'),
        )
        self.assertEqual('(A,B);', make_tree_drawable('(A,B);'))

    @unittest.skipUnless(
        importlib.util.find_spec('playwright') is not None,
        'Playwright is not installed',
    )
    @unittest.skipUnless(shutil.which('shiptv'), 'shiptv CLI is not available')
    def test_patched_html_renders_tree_and_grid_in_browser(self):
        from playwright.sync_api import sync_playwright

        html_path = self.generate_shiptv_html("(A:0.00000,B:0.00000):0.00000;")
        patch_shiptv_html_file(html_path)
        page_errors = []
        console_errors = []

        with sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            page = browser.new_page(viewport={'width': 1280, 'height': 900})
            page.on('pageerror', lambda error: page_errors.append(str(error)))
            page.on('console', lambda message: console_errors.append(message.text) if message.type == 'error' else None)
            page.goto(html_path.as_uri())
            page.wait_for_selector('#tree-plot canvas')
            page.wait_for_selector('#metadata-grid .ag-row')

            canvas_box = page.locator('#tree-plot canvas').bounding_box()
            self.assertGreater(canvas_box['width'], 0)
            self.assertGreater(canvas_box['height'], 0)
            self.assertGreaterEqual(page.locator('#metadata-grid .ag-row').count(), 2)

            page.fill('#findLeavesInput', 'A')
            page.wait_for_timeout(250)
            self.assertGreaterEqual(page.locator('#metadata-grid .ag-row-selected').count(), 1)

            page.fill('#findLeavesInput', 'B')
            page.wait_for_timeout(250)
            self.assertGreaterEqual(page.locator('#metadata-grid .ag-row-selected').count(), 1)

            page.set_viewport_size({'width': 640, 'height': 700})
            page.wait_for_timeout(250)
            resized_box = page.locator('#tree-plot canvas').bounding_box()
            self.assertGreater(resized_box['width'], 0)
            self.assertGreater(resized_box['height'], 0)
            browser.close()

        self.assertEqual([], page_errors)
        self.assertNotIn('agGrid.Grid is not a constructor', '\n'.join(console_errors))

    @unittest.skipUnless(shutil.which('shiptv'), 'shiptv CLI is not available')
    def test_non_zero_branch_lengths_are_loaded_without_normalizing_source_newick(self):
        html_path = self.generate_shiptv_html("(A:0.00000,B:0.50000):0.00000;")
        patch_shiptv_html_file(html_path)

        html = html_path.read_text(encoding='utf-8')

        self.assertIn('var newick_string = "(A:0.00000,B:0.50000):0.00000;";', html)
        self.assertIn('tree.load(makeTreeDrawable(newick_string));', html)

    @staticmethod
    def extract_shiptv_app_code(html):
        start = html.index('  var DEFAULT_GRID_COL_OPTS = {')
        end = html.rindex('</script>')
        return html[start:end]

    @staticmethod
    def generate_shiptv_html(newick):
        temp_dir = Path(tempfile.mkdtemp())
        newick_path = temp_dir / 'tree.nwk'
        metadata_path = temp_dir / 'metadata.tsv'
        html_path = temp_dir / 'tree.html'
        newick_path.write_text(newick, encoding='utf-8')
        metadata_path.write_text(
            'genome\tgroup\tscore\tcategory\n'
            'A\talpha\t1\tleft\n'
            'B\tbeta\t2\tright\n',
            encoding='utf-8',
        )
        subprocess.run(
            [
                'shiptv',
                '--newick', str(newick_path),
                '--metadata', str(metadata_path),
                '--output-html', str(html_path),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return html_path


class ExternalResourceParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.resources = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'script' and attrs.get('src'):
            self.resources.append(attrs['src'])
        if tag == 'link' and attrs.get('href'):
            self.resources.append(attrs['href'])

    @classmethod
    def find_resources(cls, html):
        parser = cls()
        parser.feed(html)
        return parser.resources


class GenbankDownloadCommandTests(SimpleTestCase):
    def test_build_download_genbank_command_passes_paths_as_arguments(self):
        command = build_download_genbank_command(
            'https://example.test/file with spaces.gbk.gz',
            '/tmp/output file.gbk',
        )

        self.assertEqual('bash', command[0])
        self.assertEqual('-o', command[1])
        self.assertEqual('pipefail', command[2])
        self.assertEqual('wget -qO- "$1" | gzip -d > "$2"', command[4])
        self.assertEqual('https://example.test/file with spaces.gbk.gz', command[6])
        self.assertEqual('/tmp/output file.gbk', command[7])

    def test_run_genbank_download_command_delegates_to_runner_with_check_enabled(self):
        class Result:
            returncode = 0

        command = ['bash', '-c', 'true']
        with patch(
            'external_tools.genbank_download_clinker_synteny.run_external_command',
            return_value=Result(),
        ) as run_command:
            returncode = run_genbank_download_command(command, timeout=25)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(command, timeout=25, shell=False, check=True)
