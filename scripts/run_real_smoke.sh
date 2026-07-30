#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
backend="${1:-}"
case "${backend}" in
  esm-if1)
    prefix="${repo_root}/.tools/envs/esm_if1"
    script="${repo_root}/scripts/smoke_esm_if1.py"
    ;;
  proteinmpnn)
    prefix="${repo_root}/.tools/envs/ligandmpnn"
    script="${repo_root}/scripts/smoke_proteinmpnn.py"
    ;;
  ligandmpnn)
    prefix="${repo_root}/.tools/envs/ligandmpnn"
    script="${repo_root}/scripts/smoke_ligandmpnn.py"
    ;;
  *)
    echo "Usage: $0 {esm-if1|proteinmpnn|ligandmpnn}" >&2
    exit 2
    ;;
esac

if [[ ! -x "${prefix}/bin/python" ]]; then
  echo "NOT_RUN: missing isolated environment ${prefix}" >&2
  exit 20
fi
conda run -p "${prefix}" python "${script}"

