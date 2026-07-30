#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
analysis_prefix="${repo_root}/.tools/envs/analysis"
if [[ ! -x "${analysis_prefix}/bin/python" ]]; then
  echo "ERROR: run make bootstrap first" >&2
  exit 2
fi
conda run -p "${analysis_prefix}" cas13-if benchmark \
  --config "${repo_root}/configs/benchmark_experimental.yaml"

