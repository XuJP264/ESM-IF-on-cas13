#!/usr/bin/env python
"""Atomically record a completed Atlas asset in its source manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path


def _sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def finalize(manifest: Path, asset: Path) -> dict[str, object]:
    if not manifest.is_file():
        raise FileNotFoundError(f"Atlas manifest is missing: {manifest}")
    if not asset.is_file():
        raise FileNotFoundError(f"completed Atlas asset is missing: {asset}")
    original = manifest.read_text(encoding="utf-8")
    match = re.search(r"^content_length:\s*(\d+)\s*$", original, flags=re.MULTILINE)
    if match is None:
        raise ValueError("Atlas manifest has no numeric content_length")
    expected_size = int(match.group(1))
    actual_size = asset.stat().st_size
    if actual_size != expected_size:
        raise ValueError(
            f"Atlas size mismatch: {actual_size} != expected {expected_size}"
        )
    digest = _sha256(asset)
    downloaded_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    replacements = {
        "size_bytes": str(actual_size),
        "sha256": digest,
        "downloaded_at": downloaded_at,
        "status": "downloaded_verified",
    }
    updated = original
    for key, value in replacements.items():
        updated, count = re.subn(
            rf"^{re.escape(key)}:.*$",
            f"{key}: {value}",
            updated,
            count=1,
            flags=re.MULTILINE,
        )
        if count != 1:
            raise ValueError(f"Atlas manifest field occurs {count} times: {key}")
    temporary = manifest.with_name(f".{manifest.name}.part-{os.getpid()}")
    temporary.write_text(updated, encoding="utf-8")
    temporary.replace(manifest)
    result: dict[str, object] = {
        "asset": str(asset),
        "manifest": str(manifest),
        "size_bytes": actual_size,
        "sha256": digest,
        "downloaded_at": downloaded_at,
        "status": "downloaded_verified",
        "is_mock": False,
    }
    return result


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=repo / "data/manifests/atlas_v1.0.yaml",
    )
    parser.add_argument(
        "--asset",
        type=Path,
        default=repo / "data/raw/atlas/v1.0/crispr-cas-atlas-v1.0.json",
    )
    arguments = parser.parse_args()
    result = finalize(arguments.manifest.resolve(), arguments.asset.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
