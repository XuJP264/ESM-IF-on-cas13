#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bundle="$(find "${repo_root}/artifacts/bundles" -mindepth 1 -maxdepth 1 -type d -name 'gpu-bundle-*' | sort | tail -n 1)"
if [[ -z "${bundle}" ]]; then
  echo "ERROR: no GPU bundle; run scripts/export_gpu_bundle.sh" >&2
  exit 2
fi
required=(bundle-manifest.json SHA256SUMS model-manifest.yaml third_party-manifest.yaml schemas/output_schema.json)
for relative in "${required[@]}"; do
  [[ -f "${bundle}/${relative}" ]] || { echo "ERROR: missing ${relative}" >&2; exit 3; }
done
(cd "${bundle}" && sha256sum --check SHA256SUMS)
python3 -m json.tool "${bundle}/bundle-manifest.json" >/dev/null
echo "VERIFIED ${bundle}"
