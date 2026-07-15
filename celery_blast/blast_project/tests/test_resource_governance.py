from pathlib import Path

from django.conf import settings
from django.core.exceptions import ValidationError
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import SimpleTestCase, override_settings

from blast_project.forms import BlastSettingsFormForward
from blast_project.models import BlastSettings
from celery_blast.resource_governance import capped_blast_threads, validate_uploaded_file_size, workflow_job_cores
from external_tools.forms import EntrezSearchForm
from external_tools.tasks import build_rpsblast_command


class ResourceGovernanceTests(SimpleTestCase):
    def test_thread_values_above_maximum_are_rejected_by_form(self):
        form = BlastSettingsFormForward(data={
            'fw_e_value': 0.001,
            'fw_word_size': 3,
            'fw_num_alignments': 100,
            'fw_max_target_seqs': 100,
            'fw_num_threads': settings.CATHI_EFFECTIVE_BLAST_THREADS + 1,
            'fw_max_hsps': 100,
        })

        self.assertFalse(form.is_valid())
        self.assertIn('fw_num_threads', form.errors)

    def test_entrez_record_count_above_maximum_is_rejected(self):
        form = EntrezSearchForm(data={
            'entrez_query': 'ribosomes AND cyanobacteria',
            'database': 'protein',
            'number_records': settings.CATHI_MAX_ENTREZ_RECORDS + 1,
        })

        self.assertFalse(form.is_valid())
        self.assertIn('number_records', form.errors)

    @override_settings(CATHI_MAX_UPLOAD_BYTES=4)
    def test_uploaded_file_above_maximum_is_rejected(self):
        uploaded_file = SimpleUploadedFile('query.faa', b'>x\nAAAA\n')

        with self.assertRaises(ValidationError):
            validate_uploaded_file_size(uploaded_file)

    def test_saved_blast_settings_are_capped_when_serialized(self):
        blast_settings = BlastSettings(
            e_value=0.001,
            word_size=3,
            num_threads=99,
            num_alignments=settings.CATHI_MAX_NUM_ALIGNMENTS + 1,
            max_target_seqs=settings.CATHI_MAX_TARGET_SEQS + 1,
            max_hsps=settings.CATHI_MAX_HSPS + 1,
        )

        values = blast_settings.values_as_fw_or_bw_dict('fw')

        self.assertEqual(str(settings.CATHI_EFFECTIVE_BLAST_THREADS), values['fw_num_threads'])
        self.assertEqual(str(settings.CATHI_MAX_NUM_ALIGNMENTS), values['fw_num_alignments'])
        self.assertEqual(str(settings.CATHI_MAX_TARGET_SEQS), values['fw_max_target_seqs'])
        self.assertEqual(str(settings.CATHI_MAX_HSPS), values['fw_max_hsps'])

    @override_settings(CATHI_WORKFLOW_JOB_CORES=2, CATHI_MAX_BLAST_THREADS=1)
    def test_effective_threads_never_exceed_job_budget_or_thread_limit(self):
        self.assertEqual(2, workflow_job_cores())
        self.assertEqual(1, capped_blast_threads(8))

    def test_rpsblast_command_caps_threads(self):
        command = build_rpsblast_command(
            'query.faa',
            'cdd',
            'out.tsv',
            {
                'rps_e_value': 0.001,
                'rps_num_threads': 99,
                'rps_max_hsps': 10,
                'rps_num_alignments': 10,
            },
        )

        self.assertEqual(str(settings.CATHI_EFFECTIVE_BLAST_THREADS), command[command.index('-num_threads') + 1])

    def test_short_and_long_tasks_route_to_different_queues(self):
        routes = settings.CELERY_TASK_ROUTES

        self.assertEqual(
            settings.CELERY_WORKFLOW_QUEUE,
            routes['blast_project.tasks.execute_reciprocal_blast_project']['queue'],
        )
        self.assertEqual(
            settings.CELERY_INTERACTIVE_QUEUE,
            routes['external_tools.tasks.entrez_search_task']['queue'],
        )
        self.assertEqual(
            settings.CELERY_MAINTENANCE_QUEUE,
            routes['external_tools.tasks.setup_cathi_download_cdd_refseq_genbank_assembly_files']['queue'],
        )

    def test_compose_workers_have_explicit_concurrency_and_runtime_limits(self):
        for compose_file in ['docker-compose.yml', 'docker-compose-production.yml']:
            candidates = [
                Path(settings.BASE_DIR).joinpath(compose_file),
                Path(settings.BASE_DIR).parent.joinpath(compose_file),
                Path.cwd().joinpath(compose_file),
            ]
            compose_path = next((path for path in candidates if path.exists()), None)
            if compose_path is None:
                self.skipTest('Compose files are not mounted in this test environment')
            content = compose_path.read_text(encoding='utf-8')

            self.assertIn('--concurrency=$${CELERY_WORKFLOW_WORKER_CONCURRENCY:-1}', content)
            self.assertIn('-Q $${CELERY_WORKFLOW_QUEUE:-long_bio}', content)
            self.assertIn('celery_worker_interactive:', content)
            self.assertIn('celery_worker_maintenance:', content)
            self.assertIn('cpus:', content)
            self.assertIn('mem_limit:', content)
            self.assertIn('pids_limit:', content)
