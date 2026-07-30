#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="${1:-all}"
export CONDA_PKGS_DIRS="${repo_root}/.tools/pkgs"
export XDG_CACHE_HOME="${repo_root}/.tools/cache"
mkdir -p "${CONDA_PKGS_DIRS}" "${XDG_CACHE_HOME}"
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is required on the GPU node" >&2
  exit 2
fi
nvidia-smi

create_env() {
  local name="$1"
  local spec="$2"
  local prefix="${repo_root}/.tools/envs/${name}"
  if [[ -x "${prefix}/bin/python" ]]; then
    conda env update -p "${prefix}" -f "${repo_root}/${spec}" --prune
  else
    conda env create -p "${prefix}" -f "${repo_root}/${spec}"
  fi
}

case "${target}" in
  all)
    create_env esm_if1 envs/esm_if1.yml
    create_env ligandmpnn envs/ligandmpnn.yml
    ;;
  esm-if1) create_env esm_if1 envs/esm_if1.yml ;;
  ligandmpnn) create_env ligandmpnn envs/ligandmpnn.yml ;;
  *) echo "ERROR: target must be all, esm-if1, or ligandmpnn" >&2; exit 3 ;;
esac
