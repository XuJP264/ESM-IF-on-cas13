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
fetch_open_pdf zhang2018_cas13d \
  https://escholarship.org/content/qt2n23w39j/qt2n23w39j.pdf
fetch_open_pdf dauparas2025_ligandmpnn \
  https://www.bakerlab.org/wp-content/uploads/2025/03/s41592-025-02626-1.pdf
fetch_open_pdf kamisetty2013_gremlin \
  https://www.bakerlab.org/wp-content/uploads/2015/12/Kamisetty_PNAS_2013.pdf

echo "Metadata-only/paywalled entries were not downloaded; see references/manifest.yaml."

