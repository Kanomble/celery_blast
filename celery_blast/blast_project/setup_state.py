from celery_blast.resource_governance import ACTIVE_TASK_STATUSES
from external_tools.models import DomainDatabase


def get_cathi_setup_state():
    domain_databases = list(DomainDatabase.objects.all())
    if not domain_databases:
        return {
            'domain_database': None,
            'setup_running': False,
        }

    setup_running = any(
        not domain_database.domain_database_loaded
        and domain_database.domain_database_download_task_result is not None
        and domain_database.domain_database_download_task_result.status in ACTIVE_TASK_STATUSES
        for domain_database in domain_databases
    )
    return {
        'domain_database': domain_databases[0],
        'setup_running': setup_running,
    }


def cathi_setup_is_running():
    return get_cathi_setup_state()['setup_running']
