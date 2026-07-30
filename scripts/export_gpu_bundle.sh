#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bundle_root="${repo_root}/artifacts/bundles"
mkdir -p "${bundle_root}"
commit="$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf 'uncommitted')"
short_commit="$(printf '%s' "${commit}" | cut -c1-12)"
bundle="${bundle_root}/gpu-bundle-${short_commit}"
if [[ -e "${bundle}" ]]; then
  echo "ERROR: refusing to overwrite ${bundle}" >&2
  exit 3
fi
mkdir -p "${bundle}/configs" "${bundle}/envs" "${bundle}/schemas"
cp "${repo_root}"/configs/*.example.yaml "${bundle}/configs/"
cp "${repo_root}"/configs/benchmark_experimental.yaml "${bundle}/configs/"
cp "${repo_root}"/envs/*.yml "${bundle}/envs/"
cp "${repo_root}"/envs/locks/* "${bundle}/envs/" 2>/dev/null || true
cp "${repo_root}/third_party/manifest.yaml" "${bundle}/third_party-manifest.yaml"
cp "${repo_root}/models/manifest.yaml" "${bundle}/model-manifest.yaml"
cp "${repo_root}/data/manifests/README.md" "${bundle}/data-manifest-readme.md"
cp "${repo_root}/src/cas13_if/refold/output_schema.json" "${bundle}/schemas/"

python3 - "${repo_root}" "${bundle}" "${commit}" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
bundle = Path(sys.argv[2])
commit = sys.argv[3]
assets = []
for relative in ("models/manifest.yaml", "third_party/manifest.yaml", "references/manifest.yaml"):
    path = root / relative
    assets.append({"path": relative, "present": path.is_file(), "size_bytes": path.stat().st_size if path.is_file() else None})
missing = []
for relative in (
    "models/esm_if1/esm_if1_gvp4_t16_142M_UR50.pt",
    "data/raw/atlas/v1.0/crispr-cas-atlas-v1.0.json",
    "data/experimental_structures/6e9f.cif",
    "data/experimental_structures/5xwp.cif",
):
    path = root / relative
    if not path.is_file():
        missing.append(relative)
manifest = {
    "schema_version": "1.0",
    "git_commit": commit,
    "is_mock": False,
    "large_assets_embedded": False,
    "manifests": assets,
    "missing_assets": missing,
    "transfer": "Use scripts/sync_assets.sh or checksum-aware rsync/rclone.",
}
(bundle / "bundle-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
PY
(cd "${bundle}" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
echo "${bundle}"
