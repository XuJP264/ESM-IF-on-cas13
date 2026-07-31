#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
paper_dir="${repo_root}/references/papers"
mkdir -p "${paper_dir}"

fetch_open_pdf() {
  local id="$1"
  local url="$2"
  local destination="${paper_dir}/${id}.pdf"
  local temporary="${destination}.part"
  if [[ -s "${destination}" ]]; then
    echo "EXISTS ${destination}"
    sha256sum "${destination}"
    return
  fi
  echo "FETCH ${url} -> ${destination}"
  curl --fail --location --retry 3 --continue-at - --output "${temporary}" "${url}"
  if [[ "$(head -c 4 "${temporary}")" != "%PDF" ]]; then
    echo "ERROR: downloaded content is not a PDF: ${url}" >&2
    exit 3
  fi
  mv "${temporary}" "${destination}"
  sha256sum "${destination}"
}

fetch_open_pdf hsu2022_esm_if1 \
  https://proceedings.mlr.press/v162/hsu22a/hsu22a.pdf
fetch_open_pdf dauparas2025_ligandmpnn \
  https://www.ipd.uw.edu/publication-pdfs/331/b896bbdf83798df6853c60bf2f2a0928/s41592-025-02626-1-3.pdf

echo "Other entries remain metadata-only: their licenses or stable direct-PDF endpoints"
echo "were not sufficiently automatable; see references/manifest.yaml."
