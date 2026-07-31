#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 SOURCE_ROOT TARGET_ROOT BUNDLE_DIR" >&2
  exit 2
fi
source_root="$1"
target_root="$2"
bundle_dir="$3"
if [[ ! -f "${bundle_dir}/ASSET_SHA256SUMS" ]]; then
  echo "ERROR: missing ${bundle_dir}/ASSET_SHA256SUMS" >&2
  exit 3
fi
for directory in models data/raw data/experimental_structures; do
  mkdir -p "${target_root}/${directory}"
  rsync -a --partial --info=progress2 "${source_root}/${directory}/" "${target_root}/${directory}/"
done
(
  cd "${target_root}"
  sha256sum --check "${bundle_dir}/ASSET_SHA256SUMS"
)
echo "ASSETS VERIFIED ${target_root}"
