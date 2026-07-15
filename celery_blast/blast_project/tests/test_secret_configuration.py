import os
import tempfile
from unittest.mock import patch

from django.core.exceptions import ImproperlyConfigured
from django.test import SimpleTestCase

from celery_blast import settings as cathi_settings


VALID_PRODUCTION_SETTINGS = {
    "cathi_env": "production",
    "debug": False,
    "secret_key": "runtime-secret-key-with-more-than-thirty-two-characters",
    "sql_engine": "django.db.backends.postgresql_psycopg2",
    "sql_database": "cathi",
    "sql_user": "cathi_app",
    "sql_password": "runtime-database-password",
    "sql_host": "postgres",
}


class SecretConfigurationTests(SimpleTestCase):
    def test_production_rejects_missing_secret_key(self):
        settings = VALID_PRODUCTION_SETTINGS.copy()
        settings["secret_key"] = ""

        with self.assertRaisesMessage(ImproperlyConfigured, "SECRET_KEY"):
            cathi_settings._validate_production_configuration(**settings)

    def test_production_rejects_placeholder_secret_key(self):
        settings = VALID_PRODUCTION_SETTINGS.copy()
        settings["secret_key"] = "CHANGE-ME-INVALID-PROD-SECRET-KEY"

        with self.assertRaisesMessage(ImproperlyConfigured, "SECRET_KEY"):
            cathi_settings._validate_production_configuration(**settings)

    def test_production_rejects_empty_database_password(self):
        settings = VALID_PRODUCTION_SETTINGS.copy()
        settings["sql_password"] = ""

        with self.assertRaisesMessage(ImproperlyConfigured, "SQL_PASSWORD"):
            cathi_settings._validate_production_configuration(**settings)

    def test_production_rejects_development_defaults(self):
        settings = VALID_PRODUCTION_SETTINGS.copy()
        settings.update(
            {
                "debug": True,
                "sql_engine": "django.db.backends.sqlite3",
                "sql_database": os.path.join(cathi_settings.BASE_DIR, "db.sqlite3"),
                "sql_user": "user",
                "sql_password": "password",
                "sql_host": "localhost",
            }
        )

        with self.assertRaises(ImproperlyConfigured) as context:
            cathi_settings._validate_production_configuration(**settings)

        message = str(context.exception)
        self.assertIn("DEBUG", message)
        self.assertIn("SQL_ENGINE", message)
        self.assertIn("SQL_PASSWORD", message)
        self.assertIn("SQL_HOST", message)

    def test_development_accepts_explicit_local_environment_secret(self):
        with patch.dict(os.environ, {"CATHI_TEST_SECRET": "local-dev-secret"}, clear=False):
            value = cathi_settings._read_secret_setting("CATHI_TEST_SECRET", required=True)

        self.assertEqual("local-dev-secret", value)

    def test_secret_file_takes_precedence_over_environment_value(self):
        secret_file = tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8")
        try:
            secret_file.write("file-backed-secret\n")
            secret_file.close()
            with patch.dict(
                os.environ,
                {
                    "CATHI_FILE_SECRET": "environment-secret",
                    "CATHI_FILE_SECRET_FILE": secret_file.name,
                },
                clear=False,
            ):
                value = cathi_settings._read_secret_setting("CATHI_FILE_SECRET", required=True)
        finally:
            os.unlink(secret_file.name)

        self.assertEqual("file-backed-secret", value)
