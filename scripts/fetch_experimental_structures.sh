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

download_uniprot_fasta() {
  local accession="$1"
  local destination="${destination_dir}/${accession}.fasta"
  local temporary="${destination}.part"
  if [[ ! -s "${destination}" ]]; then
    curl --fail --location --retry 3 \
      --output "${temporary}" \
      "https://rest.uniprot.org/uniprotkb/${accession}.fasta"
    if [[ ! -s "${temporary}" ]] || ! head -n 1 "${temporary}" | grep -q '^>'; then
      rm -f "${temporary}"
      printf 'UniProt FASTA is empty or malformed for %s\n' "${accession}" >&2
      return 1
    fi
    mv "${temporary}" "${destination}"
  fi
  sha256sum "${destination}"
}

# Stage 0003A expands the original Cas13a/Cas13d benchmark to four independent
# Cas13d parents and nine Cas13d scaffold-state units.  Keep the original
# Cas13a structures because they remain regression benchmarks.
for pdb_id in \
  6E9F 5XWP 6E9E 5XWY \
  6IV9 \
  9M38 9M30 9M33 9M34 \
  9M31 9M8Q; do
  download_structure "${pdb_id}" cif
  download_structure "${pdb_id}" pdb
done

download_metadata 6E9F 3
download_metadata 5XWP 3
download_metadata 6E9E 2
download_metadata 5XWY 2
download_metadata 6IV9 2
download_metadata 9M38 1
download_metadata 9M30 2
download_metadata 9M33 3
download_metadata 9M34 3
download_metadata 9M31 2
download_metadata 9M8Q 3

# B0MS50 remains available and is the unmutated full-length EsCas13d parent.
# The historical UrCas13d cross-reference A0A1C5SD84 currently returns an
# empty 200 response and is therefore not treated as a downloaded sequence.
download_uniprot_fasta B0MS50
