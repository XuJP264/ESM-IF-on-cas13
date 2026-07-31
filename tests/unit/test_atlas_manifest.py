import hashlib
import subprocess
import sys
from pathlib import Path


def test_finalize_atlas_manifest_is_size_gated_and_atomic(tmp_path: Path) -> None:
    asset = tmp_path / "atlas.json"
    asset.write_bytes(b"test")
    manifest = tmp_path / "atlas.yaml"
    manifest.write_text(
        "\n".join(
            (
                "name: fixture",
                "content_length: 4",
                "size_bytes: null",
                "sha256: null",
                "downloaded_at: null",
                "status: pending",
                "",
            )
        ),
        encoding="utf-8",
    )
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/finalize_atlas_manifest.py",
            "--manifest",
            str(manifest),
            "--asset",
            str(asset),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0
    updated = manifest.read_text(encoding="utf-8")
    assert "size_bytes: 4" in updated
    assert f"sha256: {hashlib.sha256(b'test').hexdigest()}" in updated
    assert "status: downloaded_verified" in updated
    assert not list(tmp_path.glob(".atlas.yaml.part-*"))

    asset.write_bytes(b"wrong-size")
    failed = subprocess.run(
        [
            sys.executable,
            "scripts/finalize_atlas_manifest.py",
            "--manifest",
            str(manifest),
            "--asset",
            str(asset),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert failed.returncode != 0
    assert "size mismatch" in failed.stderr
