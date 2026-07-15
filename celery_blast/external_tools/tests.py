import importlib.util
import shutil
import subprocess
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path
from unittest.mock import ANY, mock_open, patch

from django.contrib.auth.models import User
from django.contrib.auth.models import AnonymousUser
from django.conf import settings
from django.urls import reverse
from django.template.loader import render_to_string
from django.test import RequestFactory, SimpleTestCase, TestCase
from django_celery_results.models import TaskResult

from blast_project.models import BlastProject, BlastSettings
from external_tools.shiptv_html import make_tree_drawable, patch_shiptv_html_file
from external_tools.models import EntrezSearch, ExternalTools, QuerySequences
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
from refseq_transactions.models import BlastDatabase


def create_blast_settings():
    return BlastSettings.objects.create(
        e_value=0.001,
        word_size=3,
        num_alignments=10,
        max_target_seqs=10,
        num_threads=1,
        max_hsps=10,
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


class ExternalToolObjectAuthorizationTests(TestCase):
    def setUp(self):
        self.owner = User.objects.create_user('external-owner', password='test')
        self.other_user = User.objects.create_user('external-other', password='test')
        self.database = BlastDatabase.objects.create(
            database_name='external auth database',
            database_description='test database',
            assembly_entries=1,
            path_to_database_file='testfiles/databases/external-auth',
        )
        self.project = BlastProject.objects.create(
            project_title='external auth project',
            search_strategy='blastp',
            project_query_sequences='query.faa',
            project_user=self.owner,
            project_forward_settings=create_blast_settings(),
            project_backward_settings=create_blast_settings(),
            project_forward_database=self.database,
            project_backward_database=self.database,
            species_name_for_backward_blast='Species',
        )
        self.external_tools = ExternalTools.objects.create(
            associated_project=self.project,
            remote_or_local='local',
        )
        self.msa_task = TaskResult.objects.create(task_id='owner-msa-task', status='PROGRESS')
        self.phylo_task = TaskResult.objects.create(task_id='owner-phylo-task', status='SUCCESS')
        self.search_task = TaskResult.objects.create(task_id='owner-search-task', status='PENDING')
        self.download_task = TaskResult.objects.create(task_id='owner-download-task', status='STARTED')
        self.query_sequence = QuerySequences.objects.create(
            query_accession_id='WP_AUTH',
            query_accession_information='owned query',
            external_tool_for_query_sequence=self.external_tools,
            multiple_sequence_alignment_task=self.msa_task,
            phylogenetic_tree_construction_task=self.phylo_task,
        )
        self.entrez_search = EntrezSearch.objects.create(
            database='protein',
            entrez_user=self.owner,
            entrez_query='kinase[Title]',
            file_name='media/esearch_output/owner-search.tsv',
            fasta_file_name='media/esearch_output/owner-search.faa',
            search_task_result=self.search_task,
            download_task_result=self.download_task,
        )

    def ajax_get(self, url, user):
        self.client.force_login(user)
        return self.client.get(url, HTTP_X_REQUESTED_WITH='XMLHttpRequest')

    def assert_denied_without_state(self, response):
        self.assertEqual(404, response.status_code)
        body = response.content.decode(errors='ignore')
        for leaked_value in (
            'PROGRESS',
            'SUCCESS',
            'PENDING',
            'STARTED',
            'owner-msa-task',
            'owner-phylo-task',
            'owner-search-task',
            'owner-download-task',
            'owner-search.faa',
        ):
            self.assertNotIn(leaked_value, body)

    def test_owner_can_access_query_sequence_progress(self):
        msa_response = self.ajax_get(
            reverse('ajax_call_progress_msa_task', kwargs={'query_sequence_id': self.query_sequence.id}),
            self.owner,
        )
        phylo_response = self.ajax_get(
            reverse('ajax_call_progress_phylogeny_task', kwargs={'query_sequence_id': self.query_sequence.id}),
            self.owner,
        )

        self.assertEqual({'progress': 'PROGRESS'}, msa_response.json())
        self.assertEqual({'progress': 'SUCCESS'}, phylo_response.json())

    def test_other_user_cannot_access_query_sequence_progress(self):
        msa_response = self.ajax_get(
            reverse('ajax_call_progress_msa_task', kwargs={'query_sequence_id': self.query_sequence.id}),
            self.other_user,
        )
        phylo_response = self.ajax_get(
            reverse('ajax_call_progress_phylogeny_task', kwargs={'query_sequence_id': self.query_sequence.id}),
            self.other_user,
        )

        self.assert_denied_without_state(msa_response)
        self.assert_denied_without_state(phylo_response)

    def test_owner_can_access_entrez_progress(self):
        download_response = self.ajax_get(
            reverse('ajax_call_progress_entrezsearch_to_fasta', kwargs={'search_id': self.entrez_search.id}),
            self.owner,
        )
        search_response = self.ajax_get(
            reverse('ajax_entrez_search_progress', kwargs={'search_id': self.entrez_search.id}),
            self.owner,
        )

        self.assertEqual({'progress': 'STARTED'}, download_response.json())
        self.assertEqual({'data': 'PENDING'}, search_response.json())

    def test_other_user_cannot_access_entrez_progress_or_details(self):
        download_response = self.ajax_get(
            reverse('ajax_call_progress_entrezsearch_to_fasta', kwargs={'search_id': self.entrez_search.id}),
            self.other_user,
        )
        search_progress_response = self.ajax_get(
            reverse('ajax_entrez_search_progress', kwargs={'search_id': self.entrez_search.id}),
            self.other_user,
        )
        details_response = self.client.get(
            reverse('search_details', kwargs={'search_id': self.entrez_search.id})
        )

        self.assert_denied_without_state(download_response)
        self.assert_denied_without_state(search_progress_response)
        self.assert_denied_without_state(details_response)

    def test_anonymous_progress_requests_are_rejected(self):
        response = self.client.get(
            reverse('ajax_call_progress_msa_task', kwargs={'query_sequence_id': self.query_sequence.id}),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )

        self.assertEqual(302, response.status_code)
        self.assertIn('/blast_project/login', response['Location'])

    def test_neighboring_project_authorized_dashboard_still_works(self):
        self.client.force_login(self.owner)

        response = self.client.get(
            reverse(
                'synteny_dashboard',
                kwargs={'project_id': self.project.id, 'remote_or_local': 'local'},
            )
        )

        self.assertEqual(200, response.status_code)

    def test_owner_can_retrieve_generated_visualization_with_security_headers(self):
        self.client.force_login(self.owner)
        url = reverse(
            'load_synteny',
            kwargs={
                'project_id': self.project.id,
                'remote_or_local': 'local',
                'query_sequence_id': 'WP_AUTH',
            },
        )

        with patch('external_tools.views.get_html_results', return_value=['<script>renderPlot()</script>']):
            response = self.client.get(url)

        self.assertEqual(200, response.status_code)
        self.assertContains(response, '<script>renderPlot()</script>', html=False)
        self.assertEqual('nosniff', response['X-Content-Type-Options'])
        self.assertIn('sandbox allow-scripts', response['Content-Security-Policy'])
        self.assertNotIn('allow-same-origin', response['Content-Security-Policy'])

    def test_other_user_cannot_retrieve_generated_visualization(self):
        self.client.force_login(self.other_user)
        url = reverse(
            'load_synteny',
            kwargs={
                'project_id': self.project.id,
                'remote_or_local': 'local',
                'query_sequence_id': 'WP_AUTH',
            },
        )

        with patch('external_tools.views.get_html_results', return_value=['owned generated html']) as get_html:
            response = self.client.get(url)

        self.assertEqual(404, response.status_code)
        self.assertNotIn('owned generated html', response.content.decode(errors='ignore'))
        get_html.assert_not_called()


class EntrezSearchHtmlRenderingTests(SimpleTestCase):
    def test_external_metadata_html_like_values_render_as_text(self):
        with tempfile.TemporaryDirectory() as tempdir:
            output_path = Path(tempdir) / 'protein.tsv'
            output_path.write_text(
                '12"><script>alert(2)</script>\t'
                'Caption <b>bold</b>\t'
                'Title <script>alert(1)</script>\t'
                'Organism <img src=x onerror=alert(1)>\n',
                encoding='utf-8',
            )
            search = EntrezSearch(database='protein', file_name=str(output_path))

            html = search.get_paper_content()

        self.assertNotIn('<script', html.lower())
        self.assertNotIn('<img src=x onerror=alert(1)>', html)
        self.assertIn('&lt;script&gt;alert(1)&lt;/script&gt;', html)
        self.assertIn('Caption &lt;b&gt;bold&lt;/b&gt;', html)
        self.assertIn('Organism &lt;img src=x onerror=alert(1)&gt;', html)
        self.assertIn('https://www.ncbi.nlm.nih.gov/protein/12%22%3E%3Cscript%3Ealert%282%29%3C%2Fscript%3E', html)


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
        self.assertEqual(str(settings.CATHI_EFFECTIVE_BLAST_THREADS), command[12])
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
