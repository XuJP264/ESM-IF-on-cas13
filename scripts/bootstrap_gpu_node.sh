#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="${1:-all}"
export CONDA_PKGS_DIRS="${repo_root}/.tools/pkgs"
export XDG_CACHE_HOME="${repo_root}/.tools/cache"
export PIP_CACHE_DIR="${repo_root}/.tools/cache/pip"
export PYTHONNOUSERSITE=1
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
  if [[ "${name}" == "esm_if1" ]]; then
    if [[ ! -d "${repo_root}/third_party/esm" ]]; then
      echo "ERROR: pinned third_party/esm checkout is missing" >&2
      exit 4
    fi
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
}

case "${target}" in
  all)
    create_env analysis envs/analysis.yml
    create_env esm_if1 envs/esm_if1.yml
    create_env ligandmpnn envs/ligandmpnn.yml
    create_env bioinformatics envs/bioinformatics.yml
    ;;
  analysis) create_env analysis envs/analysis.yml ;;
  esm-if1) create_env esm_if1 envs/esm_if1.yml ;;
  ligandmpnn) create_env ligandmpnn envs/ligandmpnn.yml ;;
  bioinformatics) create_env bioinformatics envs/bioinformatics.yml ;;
  *)
    echo \
      "ERROR: target must be all, analysis, esm-if1, ligandmpnn, or bioinformatics" \
      >&2
    exit 3
    ;;
esac

case "${target}" in
  all|esm-if1)
    PYTHONNOUSERSITE=1 conda run -p "${repo_root}/.tools/envs/esm_if1" \
      python -c \
      'import esm, torch; print("esm", esm.__file__); print("cuda", torch.cuda.is_available())'
    ;;
esac
case "${target}" in
  all|ligandmpnn)
    PYTHONNOUSERSITE=1 conda run -p "${repo_root}/.tools/envs/ligandmpnn" \
      python -c \
      'import prody, torch; print("prody", prody.__version__); print("cuda", torch.cuda.is_available())'
    ;;
esac
