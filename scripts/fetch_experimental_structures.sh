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

download_metadata() {
  local pdb_id="$1"
  local entity_count="$2"
  local lower_id
  lower_id="$(printf '%s' "${pdb_id}" | tr '[:upper:]' '[:lower:]')"
  local entry_destination="${destination_dir}/${lower_id}.entry.json"
  if [[ ! -s "${entry_destination}" ]]; then
    curl --fail --location --retry 3 \
      --output "${entry_destination}.part" \
      "https://data.rcsb.org/rest/v1/core/entry/${pdb_id}"
    mv "${entry_destination}.part" "${entry_destination}"
  fi
  sha256sum "${entry_destination}"
  for entity_id in $(seq 1 "${entity_count}"); do
    local entity_destination="${destination_dir}/${lower_id}.entity_${entity_id}.json"
    if [[ ! -s "${entity_destination}" ]]; then
      curl --fail --location --retry 3 \
        --output "${entity_destination}.part" \
        "https://data.rcsb.org/rest/v1/core/polymer_entity/${pdb_id}/${entity_id}"
      mv "${entity_destination}.part" "${entity_destination}"
    fi
    sha256sum "${entity_destination}"
  done
}

# The two ternary structures are primary benchmark scaffolds. The matched
# binary structures provide same-study conformational-state context.
for pdb_id in 6E9F 5XWP 6E9E 5XWY; do
  download_structure "${pdb_id}" cif
  download_structure "${pdb_id}" pdb
done

download_metadata 6E9F 3
download_metadata 5XWP 3
download_metadata 6E9E 2
download_metadata 5XWY 2
