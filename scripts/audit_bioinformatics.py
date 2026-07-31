#!/usr/bin/env python
"""Record versions of the isolated bioinformatics toolchain."""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path


def _concise_version(name: str, output: str) -> str:
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    if name == "hmmer":
        return next((line for line in lines if line.startswith("# HMMER ")), output)
    if name == "infernal":
        return next((line for line in lines if line.startswith("# INFERNAL ")), output)
    if name == "tmalign":
        for line in lines:
            match = re.search(r"TM-align \(Version ([^)]+)\)", line)
            if match:
                return f"TM-align {match.group(1)}"
    return lines[0] if lines else ""


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    binary_dir = repo / ".tools/envs/bioinformatics/bin"
    commands = {
        "mmseqs2": [binary_dir / "mmseqs", "version"],
        "mafft": [binary_dir / "mafft", "--version"],
        "hmmer": [binary_dir / "hmmsearch", "-h"],
        "infernal": [binary_dir / "cmsearch", "-h"],
        "seqkit": [binary_dir / "seqkit", "version"],
        "foldseek": [binary_dir / "foldseek", "version"],
        "tmalign": [binary_dir / "TMalign"],
        "git_lfs": [binary_dir / "git-lfs", "version"],
    }
    records: dict[str, object] = {}
    text_lines: list[str] = []
    missing: list[str] = []
    for name, raw_command in commands.items():
        executable = Path(raw_command[0])
        command = [str(value) for value in raw_command]
        if not executable.is_file():
            records[name] = {
                "available": False,
                "executable": str(executable),
                "command": command,
            }
            missing.append(name)
            text_lines.extend((f"[{name}]", "MISSING", ""))
            continue
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
        )
        output = (completed.stdout + completed.stderr).strip()
        version = _concise_version(name, output)
        records[name] = {
            "available": True,
            "executable": str(executable),
            "command": command,
            "returncode": completed.returncode,
            "version_output": version,
        }
        text_lines.extend((f"[{name}]", version, ""))
    payload = {
        "schema_version": "1.0",
        "is_mock": False,
        "environment": "envs/bioinformatics.yml",
        "missing": missing,
        "tools": records,
    }
    output_dir = repo / "artifacts/system"
    (output_dir / "bioinformatics_tools.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_dir / "bioinformatics_tools.txt").write_text(
        "\n".join(text_lines),
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
