"""Fail CI when tracked files contain reusable secret assignments."""

import re
import subprocess
import sys
from pathlib import Path


BLOCKED_ENV_FILES = {".env.dev", ".env.prod", ".env.prod.db"}
SECRET_ASSIGNMENT = re.compile(
    r"^\s*(SECRET_KEY|SQL_PASSWORD|POSTGRES_PASSWORD|RABBITMQ_PASSWORD|AMQP_PASSWORD|FLOWER_BASIC_AUTH|JUPYTER_TOKEN)=(.+?)\s*$"
)
PLACEHOLDER_MARKERS = (
    "change-me",
    "changeme",
    "placeholder",
    "not-for-production",
    "not-a-real",
    "invalid",
    "example",
)


def tracked_files():
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return [Path(path) for path in result.stdout.decode("utf-8").split("\0") if path]


def is_placeholder(value):
    normalized = value.strip().strip("'\"").lower()
    return any(marker in normalized for marker in PLACEHOLDER_MARKERS)


def main():
    failures = []
    for path in tracked_files():
        normalized_path = path.as_posix()
        if normalized_path in BLOCKED_ENV_FILES:
            failures.append(f"{normalized_path}: real environment file must not be tracked")
            continue

        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue

        for line_number, line in enumerate(text.splitlines(), start=1):
            match = SECRET_ASSIGNMENT.match(line)
            if not match:
                continue

            variable_name, value = match.groups()
            if not is_placeholder(value):
                failures.append(
                    f"{normalized_path}:{line_number}: {variable_name} uses a non-placeholder assignment"
                )

    if failures:
        print("Tracked secret scan failed:")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print("Tracked secret scan passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
