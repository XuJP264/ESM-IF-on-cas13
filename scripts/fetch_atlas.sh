#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
url="https://storage.googleapis.com/crispr-cas-atlas-xy7q13lmk9/crispr-cas-atlas-v1.0.json"
destination_dir="${repo_root}/data/raw/atlas/v1.0"
destination="${destination_dir}/crispr-cas-atlas-v1.0.json"
temporary="${destination}.part"
mkdir -p "${destination_dir}" "${repo_root}/data/manifests"

headers="$(curl --fail --location --head --silent --show-error "${url}")"
content_length="$(printf '%s\n' "${headers}" | tr -d '\r' | awk 'tolower($1)=="content-length:" {value=$2} END {print value}')"
available_bytes="$(df --output=avail -B1 "${destination_dir}" | tail -n 1 | tr -d ' ')"

echo "SOURCE ${url}"
echo "DESTINATION ${destination}"
echo "CONTENT_LENGTH ${content_length:-unknown}"
echo "AVAILABLE_BYTES ${available_bytes}"

if [[ -n "${content_length}" && "${content_length}" =~ ^[0-9]+$ ]]; then
  required_bytes=$((content_length + content_length / 10))
  if (( available_bytes < required_bytes )); then
    echo "ERROR: insufficient disk; need 10% headroom above Content-Length" >&2
    exit 6
  fi
fi

if [[ ! -s "${destination}" ]]; then
  curl --fail --location --retry 5 --continue-at - --output "${temporary}" "${url}"
  if [[ -n "${content_length}" && "${content_length}" =~ ^[0-9]+$ ]]; then
    actual_bytes="$(stat --format='%s' "${temporary}")"
    if [[ "${actual_bytes}" != "${content_length}" ]]; then
      echo "ERROR: size mismatch ${actual_bytes} != ${content_length}" >&2
      exit 7
    fi
  fi
  mv "${temporary}" "${destination}"
fi

sha256sum "${destination}"
stat --format='%n %s bytes' "${destination}"
echo "LICENSE CC BY-NC 4.0; see https://github.com/Profluent-AI/CRISPR-Cas-Atlas"

