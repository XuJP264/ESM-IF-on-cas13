#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
analysis_prefix="${repo_root}/.tools/envs/analysis"
export CONDA_PKGS_DIRS="${repo_root}/.tools/pkgs"
export XDG_CACHE_HOME="${repo_root}/.tools/cache"
mkdir -p "${repo_root}/.tools/envs" "${repo_root}/.tools/pkgs" \
  "${repo_root}/.tools/cache" "${repo_root}/envs/locks"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is unavailable. Install Miniforge or micromamba outside system Python." >&2
  exit 2
fi

if [[ -x "${analysis_prefix}/bin/python" ]]; then
  conda env update -p "${analysis_prefix}" -f "${repo_root}/envs/analysis.yml" --prune
else
  conda env create -p "${analysis_prefix}" -f "${repo_root}/envs/analysis.yml"
fi

conda run -p "${analysis_prefix}" python -m pip install --no-deps -e "${repo_root}"
conda run -p "${analysis_prefix}" conda list --explicit > "${repo_root}/envs/locks/analysis-linux-64.explicit.txt"
conda run -p "${analysis_prefix}" conda list > "${repo_root}/envs/locks/analysis-conda-list.txt"
conda run -p "${analysis_prefix}" python -m pip freeze > "${repo_root}/envs/locks/analysis-pip-freeze.txt"
conda run -p "${analysis_prefix}" cas13-if --help >/dev/null
echo "Analysis environment ready: ${analysis_prefix}"
