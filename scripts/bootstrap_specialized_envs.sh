#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="${1:-all}"
export CONDA_PKGS_DIRS="${repo_root}/.tools/pkgs"
export XDG_CACHE_HOME="${repo_root}/.tools/cache"
export PIP_CACHE_DIR="${repo_root}/.tools/cache/pip"
export PYTHONNOUSERSITE=1
mkdir -p "${repo_root}/.tools/envs" "${CONDA_PKGS_DIRS}" \
  "${XDG_CACHE_HOME}" "${repo_root}/envs/locks"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is unavailable; system Python will not be modified." >&2
  exit 2
fi

create_and_lock() {
  local name="$1"
  local spec="$2"
  local prefix="${repo_root}/.tools/envs/${name}"
  if [[ -x "${prefix}/bin/python" ]]; then
    conda env update -p "${prefix}" -f "${repo_root}/${spec}" --prune
  else
    conda env create -p "${prefix}" -f "${repo_root}/${spec}"
  fi
  if [[ "${name}" == "esm_if1" ]]; then
    conda run -p "${prefix}" python -m pip install \
      --no-deps --no-build-isolation --force-reinstall \
      "${repo_root}/third_party/esm"
  fi
  conda run -p "${prefix}" python -m pip install --no-deps -e "${repo_root}"
  conda list -p "${prefix}" --explicit \
    > "${repo_root}/envs/locks/${name}-linux-64.explicit.txt"
  conda list -p "${prefix}" \
    > "${repo_root}/envs/locks/${name}-conda-list.txt"
  conda run -p "${prefix}" python -m pip freeze \
    > "${repo_root}/envs/locks/${name}-pip-freeze.txt"
  if [[ "${name}" == "bioinformatics" ]]; then
    conda run -p "${prefix}" python \
      "${repo_root}/scripts/audit_bioinformatics.py"
  fi
  echo "Environment ready: ${name} (${prefix})"
}

case "${target}" in
  all)
    create_and_lock esm_if1 envs/esm_if1.yml
    create_and_lock ligandmpnn envs/ligandmpnn.yml
    create_and_lock bioinformatics envs/bioinformatics.yml
    ;;
  esm-if1) create_and_lock esm_if1 envs/esm_if1.yml ;;
  ligandmpnn) create_and_lock ligandmpnn envs/ligandmpnn.yml ;;
  bioinformatics) create_and_lock bioinformatics envs/bioinformatics.yml ;;
  *)
    echo "ERROR: target must be all, esm-if1, ligandmpnn, or bioinformatics" >&2
    exit 3
    ;;
esac
