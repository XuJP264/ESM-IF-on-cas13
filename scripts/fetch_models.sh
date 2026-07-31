#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "${repo_root}/models/esm_if1" "${repo_root}/models/proteinmpnn" \
  "${repo_root}/models/ligandmpnn"

download_atomic() {
  local url="$1"
  local destination="$2"
  local temporary="${destination}.part"
  if [[ ! -s "${destination}" ]]; then
    echo "FETCH ${url} -> ${destination}"
    curl --fail --location --retry 3 --continue-at - --output "${temporary}" "${url}"
    mv "${temporary}" "${destination}"
  fi
  sha256sum "${destination}"
  stat --format='%n %s bytes' "${destination}"
}

download_atomic \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt \
  "${repo_root}/models/esm_if1/esm_if1_gvp4_t16_142M_UR50.pt"

if [[ ! -d "${repo_root}/third_party/ProteinMPNN/.git" \
   || ! -d "${repo_root}/third_party/LigandMPNN/.git" ]]; then
  bash "${repo_root}/scripts/fetch_third_party.sh"
fi

copy_checkpoint() {
  local source="$1"
  local destination="$2"
  if [[ ! -s "${source}" ]]; then
    echo "ERROR: pinned upstream checkpoint missing: ${source}" >&2
    exit 5
  fi
  if [[ ! -s "${destination}" ]]; then
    cp "${source}" "${destination}.part"
    mv "${destination}.part" "${destination}"
  fi
  sha256sum "${destination}"
  stat --format='%n %s bytes' "${destination}"
}

copy_checkpoint \
  "${repo_root}/third_party/ProteinMPNN/vanilla_model_weights/v_48_020.pt" \
  "${repo_root}/models/proteinmpnn/v_48_020.pt"
# LigandMPNN deliberately does not store weights in Git. Fetch only the selected,
# pinned model variants from the official IPD distribution endpoint; do not run
# the upstream all-model downloader because that would retrieve unused weights.
download_atomic \
  https://files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_020.pt \
  "${repo_root}/models/ligandmpnn/proteinmpnn_v_48_020.pt"
download_atomic \
  https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_010_25.pt \
  "${repo_root}/models/ligandmpnn/ligandmpnn_v_32_010_25.pt"
download_atomic \
  https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_020.pt \
  "${repo_root}/models/ligandmpnn/solublempnn_v_48_020.pt"
