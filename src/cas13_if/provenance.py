"""Immutable run records and file provenance."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

from cas13_if.config import ConfigDict, config_hash


class RunExistsError(FileExistsError):
    """Raised rather than overwriting an existing run."""


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def git_metadata(repo: Path) -> dict[str, Any]:
    def run(*arguments: str) -> str:
        result = subprocess.run(
            ["git", *arguments],
            cwd=repo,
            text=True,
            capture_output=True,
            check=False,
        )
        return result.stdout.strip() if result.returncode == 0 else "unavailable"

    status = run("status", "--porcelain")
    return {
        "commit": run("rev-parse", "HEAD"),
        "short_sha": run("rev-parse", "--short", "HEAD"),
        "branch": run("branch", "--show-current"),
        "remote": run("remote", "get-url", "origin"),
        "dirty": bool(status and status != "unavailable"),
        "status_porcelain": status.splitlines() if status != "unavailable" else [],
    }


def environment_snapshot() -> str:
    distributions = sorted(
        (
            distribution.name,
            distribution.version,
        )
        for distribution in importlib.metadata.distributions()
    )
    lines = [
        f"python={sys.version.replace(chr(10), ' ')}",
        f"executable={sys.executable}",
        f"platform={platform.platform()}",
    ]
    lines.extend(f"{name}=={version}" for name, version in distributions)
    return "\n".join(lines) + "\n"


def make_run_id(
    experiment: str,
    resolved_config: ConfigDict,
    short_sha: str,
    *,
    now: datetime | None = None,
) -> str:
    timestamp = now or datetime.now().astimezone()
    safe_experiment = "".join(
        character if character.isalnum() or character in "-_" else "-"
        for character in experiment.lower()
    ).strip("-")
    sha = short_sha if short_sha != "unavailable" else "nogit"
    return f"{timestamp:%Y%m%d}-{safe_experiment}-{config_hash(resolved_config)}-{sha}"


@dataclass
class RunRecorder:
    root: Path
    experiment: str
    resolved_config: ConfigDict
    command: list[str]
    repo_root: Path
    is_mock: bool
    run_dir: Path = field(init=False)
    failures: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        git = git_metadata(self.repo_root)
        run_id = make_run_id(
            self.experiment, self.resolved_config, str(git["short_sha"])
        )
        self.run_dir = self.root / run_id
        try:
            self.run_dir.mkdir(parents=True, exist_ok=False)
        except FileExistsError as exc:
            raise RunExistsError(f"refusing to overwrite run: {self.run_dir}") from exc
        atomic_write_text(
            self.run_dir / "resolved_config.yaml",
            yaml.safe_dump(self.resolved_config, sort_keys=True),
        )
        atomic_write_text(
            self.run_dir / "command.txt",
            shlex.join(self.command) + "\n",
        )
        atomic_write_text(self.run_dir / "environment.txt", environment_snapshot())
        hardware_source = self.repo_root / "artifacts/system/hardware.json"
        hardware = (
            hardware_source.read_text(encoding="utf-8")
            if hardware_source.is_file()
            else json.dumps({"status": "not_available", "is_mock": self.is_mock})
        )
        atomic_write_text(self.run_dir / "hardware.json", hardware.rstrip() + "\n")
        atomic_write_text(
            self.run_dir / "git.json",
            json.dumps(git, indent=2, sort_keys=True) + "\n",
        )
        for filename, default in (
            ("input_manifest.json", {"files": [], "is_mock": self.is_mock}),
            ("output_manifest.json", {"files": [], "is_mock": self.is_mock}),
            ("metrics.json", {"metrics": {}, "is_mock": self.is_mock}),
        ):
            atomic_write_text(
                self.run_dir / filename,
                json.dumps(default, indent=2, sort_keys=True) + "\n",
            )
        atomic_write_text(self.run_dir / "failures.jsonl", "")
        atomic_write_text(self.run_dir / "stdout.log", "")
        atomic_write_text(self.run_dir / "stderr.log", "")

    def record_failure(self, stage: str, message: str) -> None:
        failure = {"stage": stage, "message": message, "is_mock": self.is_mock}
        self.failures.append(failure)
        content = "".join(
            json.dumps(item, sort_keys=True) + "\n" for item in self.failures
        )
        atomic_write_text(self.run_dir / "failures.jsonl", content)

    def finish(
        self,
        *,
        success: bool,
        metrics: dict[str, Any] | None = None,
        outputs: list[dict[str, Any]] | None = None,
    ) -> None:
        atomic_write_text(
            self.run_dir / "metrics.json",
            json.dumps(
                {"metrics": metrics or {}, "is_mock": self.is_mock},
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
        atomic_write_text(
            self.run_dir / "output_manifest.json",
            json.dumps(
                {"files": outputs or [], "is_mock": self.is_mock},
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
        exit_code = 0 if success else 1
        atomic_write_text(self.run_dir / "exit_code", f"{exit_code}\n")
        marker = "SUCCESS" if success else "FAILED"
        atomic_write_text(self.run_dir / marker, "\n")
