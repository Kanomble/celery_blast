import os
import subprocess
import unicodedata

import yaml


class EntrezQueryValidationError(ValueError):
    pass


def validate_entrez_query(entrez_query):
    if entrez_query is None:
        return ""

    query = str(entrez_query)
    invalid_characters = [
        character
        for character in query
        if unicodedata.category(character) in {"Cc", "Zl", "Zp"}
    ]
    if invalid_characters:
        raise EntrezQueryValidationError(
            "Entrez query contains unsupported control characters or line breaks."
        )
    return query


def load_snakemake_config(config_path):
    with open(config_path, encoding="utf-8") as config_file:
        data = yaml.safe_load(config_file) or {}

    if not isinstance(data, dict):
        raise EntrezQueryValidationError("Snakemake configuration must be a mapping.")
    return data


def write_snakemake_config(config_path, data):
    with open(config_path, "w", encoding="utf-8") as config_file:
        yaml.safe_dump(
            data,
            config_file,
            allow_unicode=True,
            default_flow_style=False,
            sort_keys=False,
        )


def update_snakemake_config_entrez_query(config_path, entrez_query):
    data = load_snakemake_config(config_path)
    data["entrez_query"] = validate_entrez_query(entrez_query)
    write_snakemake_config(config_path, data)


def build_remote_blast_command(
    *,
    search_strategy,
    database,
    outfmt,
    output_path,
    word_size,
    e_value,
    num_alignments,
    query_path,
    entrez_query="",
    max_hsps=None,
    num_threads=None,
):
    command = [
        str(search_strategy),
        "-db", str(database),
        "-outfmt", str(outfmt),
        "-out", str(output_path),
        "-word_size", str(word_size),
        "-evalue", str(e_value),
        "-num_alignments", str(num_alignments),
    ]

    if max_hsps is not None:
        command.extend(["-max_hsps", str(max_hsps)])
    if num_threads is not None:
        command.extend(["-num_threads", str(num_threads)])

    command.extend(["-query", str(query_path), "-remote"])

    validated_query = validate_entrez_query(entrez_query)
    if validated_query.strip() and validated_query.strip() != "all[organism]":
        command.extend(["-entrez_query", validated_query])

    return command


def run_remote_blast_command(command, stderr_path):
    stderr_dir = os.path.dirname(str(stderr_path))
    if stderr_dir:
        os.makedirs(stderr_dir, exist_ok=True)
    with open(stderr_path, "w", encoding="utf-8") as stderr_file:
        return subprocess.run(command, stderr=stderr_file, check=True)
