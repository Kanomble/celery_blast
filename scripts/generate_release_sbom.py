"""Generate a small CycloneDX SBOM for CATHI's principal scientific tools."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def component_from_artifact(artifact: dict[str, object]) -> dict[str, object]:
    external_references = [
        {
            "type": "distribution",
            "url": artifact["url"],
        }
    ]
    hashes = [
        {
            "alg": "SHA-256",
            "content": artifact["sha256"],
        }
    ]
    return {
        "type": "application",
        "bom-ref": f"cathi-tool:{artifact['id']}",
        "name": artifact["name"],
        "version": artifact["version"],
        "hashes": hashes,
        "externalReferences": external_references,
        "properties": [
            {"name": "cathi:artifact:type", "value": artifact["type"]},
            {"name": "cathi:artifact:id", "value": artifact["id"]},
        ],
    }


def build_sbom(manifest: dict[str, object]) -> dict[str, object]:
    artifacts = manifest["artifacts"]
    components = [
        component_from_artifact(artifact)
        for artifact in artifacts
        if artifact.get("principal_scientific_tool")
    ]

    components.extend(
        {
            "type": "container",
            "bom-ref": f"cathi-base:{base['name']}",
            "name": base["name"],
            "version": base["reference"],
            "properties": [
                {"name": "cathi:base-image:usage", "value": base["usage"]},
                {"name": "cathi:base-image:platform", "value": base["platform"]},
            ],
        }
        for base in manifest["base_images"]
    )

    return {
        "bomFormat": "CycloneDX",
        "specVersion": "1.5",
        "serialNumber": "urn:uuid:00000000-0000-0000-0000-000000000000",
        "version": 1,
        "metadata": {
            "timestamp": "1970-01-01T00:00:00+00:00",
            "component": {
                "type": "application",
                "name": "CATHI runtime image scientific tools",
                "version": "source-controlled",
            },
            "properties": [
                {
                    "name": "cathi:supported-architectures",
                    "value": ",".join(manifest["supported_architectures"]),
                }
            ],
        },
        "components": components,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-manifest",
        default="build/remote-artifacts.json",
        help="Path to CATHI remote artifact manifest.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write CycloneDX JSON SBOM.",
    )
    args = parser.parse_args()

    manifest_path = Path(args.artifact_manifest)
    output_path = Path(args.output)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(build_sbom(manifest), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
