#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
required=(
  AGENTS.md .agent/PLANS.md README.md pyproject.toml Makefile
  docs/STATUS.md docs/DECISIONS.md workflow/Snakefile
  third_party/manifest.yaml models/manifest.yaml references/manifest.yaml
  artifacts/system/hardware.json
)
for relative in "${required[@]}"; do
  if [[ ! -f "${repo_root}/${relative}" ]]; then
    echo "ERROR: missing ${relative}" >&2
    exit 10
  fi
done

while IFS= read -r -d '' script; do
  bash -n "${script}"
done < <(find "${repo_root}/scripts" -maxdepth 1 -type f -name '*.sh' -print0)
echo "Shell scripts parse successfully"

python3 - "${repo_root}" <<'PY'
import json
import sys
from pathlib import Path

import yaml

root = Path(sys.argv[1])
json.loads((root / "artifacts/system/hardware.json").read_text())
for relative in (
    "third_party/manifest.yaml",
    "models/manifest.yaml",
    "references/manifest.yaml",
):
    value = yaml.safe_load((root / relative).read_text())
    if not isinstance(value, dict):
        raise SystemExit(f"invalid mapping: {relative}")
print("Static manifests valid")
PY

if [[ -x "${repo_root}/.tools/envs/analysis/bin/python" ]]; then
  conda run -p "${repo_root}/.tools/envs/analysis" ruff check \
    "${repo_root}/src" "${repo_root}/tests" "${repo_root}/scripts"
  conda run -p "${repo_root}/.tools/envs/analysis" ruff format --check \
    "${repo_root}/src" "${repo_root}/tests" "${repo_root}/scripts"
  conda run -p "${repo_root}/.tools/envs/analysis" mypy "${repo_root}/src"
  conda run -p "${repo_root}/.tools/envs/analysis" pytest \
    -m "not real_model and not network and not slow" "${repo_root}/tests"
else
  echo "ERROR: analysis environment missing; run make bootstrap" >&2
  exit 11
fi
