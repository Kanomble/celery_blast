from subprocess import TimeoutExpired
from unittest import TestCase

from celery_blast.processes import (
    ExternalCommandError,
    ExternalCommandTimeout,
    run_external_command,
)


class FakeProcess:
    def __init__(self, pid=123, returncode=0, timeout=False):
        self.pid = pid
        self.returncode = returncode
        self.timeout = timeout
        self.killed = False
        self.wait_timeout = None

    def wait(self, timeout=None):
        self.wait_timeout = timeout
        if self.timeout:
            raise TimeoutExpired(cmd='fake', timeout=timeout)
        return self.returncode

    def kill(self):
        self.killed = True


class ExternalCommandRunnerTests(TestCase):
    def test_run_external_command_returns_result(self):
        process = FakeProcess(returncode=0)
        calls = []

        def popen_factory(command, shell=False):
            calls.append((command, shell))
            return process

        result = run_external_command(
            ['snakemake', '--cores', '1'],
            timeout=30,
            shell=False,
            popen_factory=popen_factory,
        )

        self.assertEqual(0, result.returncode)
        self.assertEqual(123, result.pid)
        self.assertEqual(30, process.wait_timeout)
        self.assertEqual([(['snakemake', '--cores', '1'], False)], calls)

    def test_run_external_command_passes_streams_when_provided(self):
        process = FakeProcess(returncode=0)
        calls = []
        stdout = object()
        stderr = object()

        def popen_factory(command, **kwargs):
            calls.append((command, kwargs))
            return process

        result = run_external_command(
            ['mafft', 'input.faa'],
            timeout=30,
            stdout=stdout,
            stderr=stderr,
            popen_factory=popen_factory,
        )

        self.assertEqual(0, result.returncode)
        self.assertEqual(
            [(['mafft', 'input.faa'], {'shell': False, 'stdout': stdout, 'stderr': stderr})],
            calls,
        )

    def test_run_external_command_raises_on_nonzero_returncode_when_checking(self):
        process = FakeProcess(returncode=2)

        with self.assertRaises(ExternalCommandError) as exc:
            run_external_command(
                ['snakemake'],
                timeout=30,
                check=True,
                popen_factory=lambda command, shell=False: process,
            )

        self.assertEqual(2, exc.exception.returncode)

    def test_run_external_command_kills_process_tree_on_timeout(self):
        process = FakeProcess(timeout=True)
        killed = []

        with self.assertRaises(ExternalCommandTimeout):
            run_external_command(
                ['snakemake'],
                timeout=30,
                popen_factory=lambda command, shell=False: process,
                process_tree_killer=lambda current_process: killed.append(current_process.pid),
            )

        self.assertEqual([123], killed)

    def test_run_external_command_kills_process_tree_on_cleanup_exception(self):
        class SoftLimit(Exception):
            pass

        class SoftLimitProcess(FakeProcess):
            def wait(self, timeout=None):
                raise SoftLimit()

        process = SoftLimitProcess()
        killed = []

        with self.assertRaises(SoftLimit):
            run_external_command(
                ['snakemake'],
                timeout=30,
                cleanup_exceptions=(SoftLimit,),
                popen_factory=lambda command, shell=False: process,
                process_tree_killer=lambda current_process: killed.append(current_process.pid),
            )

        self.assertEqual([123], killed)
