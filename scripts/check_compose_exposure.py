"""Validate Compose port exposure and administrative service settings."""

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def compose_config(*args):
    command = ["docker", "compose", *args, "config", "--no-env-resolution", "--format", "json"]
    result = subprocess.run(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(result.stderr.strip() or "docker compose config failed")
    return json.loads(result.stdout)


def service_networks(service):
    networks = service.get("networks", {})
    if isinstance(networks, dict):
        return set(networks)
    return set(networks)


def ports_for(service):
    return service.get("ports") or []


def host_ip(port):
    return port.get("host_ip") or port.get("host_ip".replace("_", ""))


def command_text(service):
    command = service.get("command") or ""
    if isinstance(command, list):
        return " ".join(str(part) for part in command)
    return str(command)


def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)


def check_production_default_exposure(config):
    services = config["services"]
    for name in ("flower", "jupyter_notebook"):
        assert_true(name not in services, f"{name} must be absent unless the admin profile is enabled")

    published = {name: ports_for(service) for name, service in services.items() if ports_for(service)}

    assert_true(set(published) == {"nginx"}, "default production must publish only nginx")
    nginx_ports = published["nginx"]
    assert_true(len(nginx_ports) == 1, "nginx must publish exactly one host port")
    assert_true(str(nginx_ports[0].get("target")) == "8080", "nginx must publish unprivileged container port 8080")
    assert_true(not ports_for(services["web"]), "web must not publish a direct host port")
    assert_true(not ports_for(services["postgres"]), "postgres must not publish a production host port")
    assert_true(not ports_for(services["rabbitmq"]), "rabbitmq must not publish a production host port")


def check_service_reachability(config):
    services = config["services"]
    for app_service in ("web", "celery_worker"):
        app_networks = service_networks(services[app_service])
        assert_true("data" in app_networks, f"{app_service} must join the data network")

    assert_true("data" in service_networks(services["postgres"]), "postgres must join the data network")
    assert_true("data" in service_networks(services["rabbitmq"]), "rabbitmq must join the data network")
    assert_true("app" in service_networks(services["web"]), "web must join the app network")
    assert_true("app" in service_networks(services["nginx"]), "nginx must reach web on the app network")
    assert_true("edge" in service_networks(services["nginx"]), "nginx must join the edge network")


def check_admin_services(config):
    services = config["services"]
    for name in ("flower", "jupyter_notebook"):
        assert_true(services[name].get("profiles") == ["admin"], f"{name} must be admin-profile only")
        admin_ports = ports_for(services[name])
        assert_true(admin_ports, f"{name} must publish a loopback admin port when enabled")
        assert_true(
            all((port.get("host_ip") or "") == "127.0.0.1" for port in admin_ports),
            f"{name} admin ports must bind to 127.0.0.1",
        )

    flower_command = command_text(services["flower"])
    assert_true("FLOWER_BASIC_AUTH" in flower_command, "flower must validate FLOWER_BASIC_AUTH")
    assert_true("--basic_auth" in flower_command, "flower must enable basic authentication")

    jupyter_command = command_text(services["jupyter_notebook"])
    assert_true("JUPYTER_TOKEN" in jupyter_command, "jupyter must validate JUPYTER_TOKEN")
    assert_true("--NotebookApp.token" in jupyter_command, "jupyter must set an explicit token")
    assert_true("--NotebookApp.password=" not in jupyter_command, "jupyter must not force a blank password")


def check_development_loopback(config):
    services = config["services"]
    for name in ("postgres", "rabbitmq", "flower"):
        service_ports = ports_for(services[name])
        assert_true(service_ports, f"development {name} must keep an explicit local port")
        assert_true(
            all((port.get("host_ip") or "") == "127.0.0.1" for port in service_ports),
            f"development {name} ports must bind to 127.0.0.1",
        )


def check_rabbitmq_consumer_timeout(config):
    rabbitmq = config["services"]["rabbitmq"]
    command = command_text(rabbitmq)
    volumes = rabbitmq.get("volumes") or []
    env_files = rabbitmq.get("env_file") or []
    env = rabbitmq.get("environment") or {}

    assert_true(
        "RABBITMQ_CONSUMER_TIMEOUT_MS" in command,
        "rabbitmq must generate consumer_timeout from RABBITMQ_CONSUMER_TIMEOUT_MS",
    )
    assert_true(
        "RABBITMQ_CONFIG_FILE=/tmp/rabbitmq" in command,
        "rabbitmq must start with the generated RabbitMQ config file",
    )
    assert_true(
        any("rabbitmq.conf.template" in str(volume) for volume in volumes),
        "rabbitmq must mount the consumer timeout config template",
    )
    assert_true(
        env_files or "RABBITMQ_CONSUMER_TIMEOUT_MS" in env,
        "rabbitmq must receive RABBITMQ_CONSUMER_TIMEOUT_MS from env config",
    )


def main():
    try:
        production = compose_config("-f", "docker-compose-production.yml")
        production_admin = compose_config("-f", "docker-compose-production.yml", "--profile", "admin")
        check_production_default_exposure(production)
        check_service_reachability(production)
        check_rabbitmq_consumer_timeout(production)
        check_admin_services(production_admin)

        development = compose_config("-f", "docker-compose.yml")
        check_development_loopback(development)
        check_rabbitmq_consumer_timeout(development)
    except AssertionError as exc:
        print(f"Compose exposure check failed: {exc}")
        return 1

    print("Compose exposure check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
