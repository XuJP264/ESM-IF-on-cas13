#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bundle="${1:-}"
asset_root="${2:-}"
if [[ -z "${bundle}" ]]; then
  latest_manifest="$(
    find "${repo_root}/artifacts/bundles" -mindepth 2 -maxdepth 2 \
      -type f -path '*/gpu-bundle-*/bundle-manifest.json' \
      -printf '%T@ %p\n' |
      sort -n |
      tail -n 1 |
      cut -d ' ' -f 2-
  )"
  bundle="${latest_manifest%/bundle-manifest.json}"
fi
if [[ -z "${bundle}" ]]; then
  echo "ERROR: no GPU bundle; run scripts/export_gpu_bundle.sh" >&2
  exit 2
fi
bundle="$(cd "${bundle}" && pwd)"
required=(
  ASSET_SHA256SUMS
  bundle-manifest.json
  repo-clone.txt
  SHA256SUMS
  model-manifest.yaml
  third_party-manifest.yaml
  schemas/output_schema.json
  scripts/bootstrap_gpu_node.sh
  scripts/launch_gpu_tmux.sh
  scripts/sync_assets.sh
)
for relative in "${required[@]}"; do
  [[ -f "${bundle}/${relative}" ]] || { echo "ERROR: missing ${relative}" >&2; exit 3; }
done
if [[ -d "${bundle}/inputs" ]] && find "${bundle}/inputs" -type f -print -quit | grep -q .; then
  stage_required=(
    inputs/manifests/all_jobs.jsonl
    inputs/manifests/candidate_inventory.csv
    inputs/manifests/summary.json
    inputs/manifests/SHA256SUMS
    inputs/expected_outputs/expected_outputs.jsonl
    inputs/retry_manifests/failed_jobs.jsonl
  )
  for relative in "${stage_required[@]}"; do
    [[ -f "${bundle}/${relative}" ]] || {
      echo "ERROR: missing Stage 0003 input ${relative}" >&2
      exit 3
    }
  done
fi
(cd "${bundle}" && sha256sum --check SHA256SUMS)
python3 - "${bundle}/bundle-manifest.json" <<'PY'
import json
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text())
required = {
    "schema_version",
    "git_commit",
    "git_worktree_dirty_at_export",
    "config_hash",
    "assets",
    "missing_assets",
    "input_files",
}
missing = sorted(required.difference(manifest))
if missing:
    raise SystemExit(f"bundle manifest fields missing: {missing}")
if manifest["large_assets_embedded"]:
    raise SystemExit("large assets must not be embedded in the GPU bundle")
print(
    "Bundle manifest valid:",
    manifest["git_commit"],
    "dirty=" + str(manifest["git_worktree_dirty_at_export"]).lower(),
)
PY
if [[ -n "${asset_root}" ]]; then
  (
    cd "${asset_root}"
    sha256sum --check "${bundle}/ASSET_SHA256SUMS"
  )
else
  echo "NOTICE: asset bytes not checked; pass ASSET_ROOT as argument 2 after sync"
fi
echo "VERIFIED ${bundle}"
