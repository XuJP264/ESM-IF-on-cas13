#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

fetch_repo() {
  local name="$1"
  local url="$2"
  local commit="$3"
  local destination="${repo_root}/third_party/${name}"
  if [[ ! -d "${destination}/.git" ]]; then
    git clone --filter=blob:none "${url}" "${destination}"
  fi
  git -C "${destination}" fetch --depth 1 origin "${commit}"
  git -C "${destination}" checkout --detach "${commit}"
  local actual
  actual="$(git -C "${destination}" rev-parse HEAD)"
  if [[ "${actual}" != "${commit}" ]]; then
    echo "ERROR: ${name} resolved to ${actual}, expected ${commit}" >&2
    exit 3
  fi
  if [[ -n "$(git -C "${destination}" status --porcelain)" ]]; then
    echo "ERROR: refusing modified third-party checkout: ${destination}" >&2
    exit 4
  fi
  echo "${name} ${actual}"
}

fetch_repo esm https://github.com/facebookresearch/esm.git \
  2b369911bb5b4b0dda914521b9475cad1656b2ac
fetch_repo SynTnpBs https://github.com/pyskop/SynTnpBs.git \
  f3ea8e69c6f71baa56c4bb388e9df0489720f968
fetch_repo ProteinMPNN https://github.com/dauparas/ProteinMPNN.git \
  8907e6671bfbfc92303b5f79c4b5e6ce47cdef57
fetch_repo LigandMPNN https://github.com/dauparas/LigandMPNN.git \
  26ec57ac976ade5379920dbd43c7f97a91cf82de

