#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
destination_dir="${repo_root}/data/experimental_structures"
mkdir -p "${destination_dir}"

download_structure() {
  local pdb_id="$1"
  local extension="$2"
  local lower_id
  lower_id="$(printf '%s' "${pdb_id}" | tr '[:upper:]' '[:lower:]')"
  local destination="${destination_dir}/${lower_id}.${extension}"
  local temporary="${destination}.part"
  local url="https://files.rcsb.org/download/${pdb_id}.${extension}"
  if [[ ! -s "${destination}" ]]; then
    curl --fail --location --retry 3 --continue-at - --output "${temporary}" "${url}"
    mv "${temporary}" "${destination}"
  fi
  sha256sum "${destination}"
}

for pdb_id in 6E9F 5XWP; do
  download_structure "${pdb_id}" cif
  download_structure "${pdb_id}" pdb
done

