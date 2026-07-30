#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 SOURCE_ROOT TARGET_ROOT" >&2
  exit 2
fi
source_root="$1"
target_root="$2"
for directory in models data/raw data/experimental_structures; do
  mkdir -p "${target_root}/${directory}"
  rsync -a --partial --info=progress2 "${source_root}/${directory}/" "${target_root}/${directory}/"
done

