from dataclasses import dataclass
from subprocess import Popen, SubprocessError, TimeoutExpired

try:
    import psutil
except ImportError:  # pragma: no cover - psutil is available in the runtime image.
    psutil = None


@dataclass(frozen=True)
class CommandResult:
    command: object
    returncode: int
    pid: int


class ExternalCommandError(Exception):
    def __init__(self, command, returncode=None, message=None, pid=None):
        self.command = command
        self.returncode = returncode
        self.pid = pid
        if message is None:
            message = "external command failed"
            if returncode is not None:
                message = "{} with return code {}".format(message, returncode)
        super().__init__(message)


class ExternalCommandTimeout(ExternalCommandError):
    pass


def kill_process_tree(process):
    if psutil is None:
        process.kill()
        return

    parent = psutil.Process(process.pid)
    for child in parent.children(recursive=True):
        child.kill()
    parent.kill()


def run_external_command(
        command,
        timeout,
        shell=False,
        logger=None,
        check=False,
        cleanup_exceptions=(),
        stdout=None,
        stderr=None,
        popen_factory=Popen,
        process_tree_killer=kill_process_tree):
    popen_kwargs = {'shell': shell}
    if stdout is not None:
        popen_kwargs['stdout'] = stdout
    if stderr is not None:
        popen_kwargs['stderr'] = stderr
    process = popen_factory(command, **popen_kwargs)
    if logger is not None:
        logger.info('waiting for popen instance {} to finish with timeout set to {}'.format(process.pid, timeout))

    try:
        returncode = process.wait(timeout=timeout)
    except TimeoutExpired as exc:
        if logger is not None:
            logger.info('timeout expired ... trying to kill process {}'.format(process.pid))
        process_tree_killer(process)
        raise ExternalCommandTimeout(
            command,
            message='external command timed out: {}'.format(command),
            pid=process.pid,
        ) from exc
    except SubprocessError as exc:
        process_tree_killer(process)
        raise ExternalCommandError(
            command,
            message='external command raised subprocess error: {}'.format(exc),
            pid=process.pid,
        ) from exc
    except cleanup_exceptions:
        process_tree_killer(process)
        raise

    if logger is not None:
        logger.info('returncode : {}'.format(returncode))

    if check and returncode != 0:
        raise ExternalCommandError(command, returncode=returncode)

    return CommandResult(command=command, returncode=returncode, pid=process.pid)
