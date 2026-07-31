#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bundle_root="${repo_root}/artifacts/bundles"
input_root="${1:-}"
mkdir -p "${bundle_root}"
commit="$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf 'uncommitted')"
short_commit="$(printf '%s' "${commit}" | cut -c1-12)"
config_hash="$(
  find "${repo_root}/configs" "${repo_root}/envs" -type f \
    \( -name '*.yaml' -o -name '*.yml' \) -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    | sha256sum \
    | cut -c1-10
)"
bundle="${bundle_root}/gpu-bundle-${short_commit}-${config_hash}"
if [[ -e "${bundle}" ]]; then
  echo "ERROR: refusing to overwrite ${bundle}" >&2
  exit 3
fi
mkdir -p \
  "${bundle}/configs" \
  "${bundle}/containers" \
  "${bundle}/docs" \
  "${bundle}/envs/locks" \
  "${bundle}/inputs" \
  "${bundle}/schemas" \
  "${bundle}/scripts"
cp "${repo_root}"/configs/*.yaml "${bundle}/configs/"
cp "${repo_root}"/containers/* "${bundle}/containers/"
cp "${repo_root}"/docs/GPU_MIGRATION.md "${bundle}/docs/"
cp "${repo_root}"/envs/*.yml "${bundle}/envs/"
cp "${repo_root}"/envs/locks/* "${bundle}/envs/locks/" 2>/dev/null || true
cp \
  "${repo_root}/scripts/bootstrap_gpu_node.sh" \
  "${repo_root}/scripts/launch_gpu_tmux.sh" \
  "${repo_root}/scripts/sync_assets.sh" \
  "${repo_root}/scripts/verify_gpu_bundle.sh" \
  "${bundle}/scripts/"
cp "${repo_root}/third_party/manifest.yaml" "${bundle}/third_party-manifest.yaml"
cp "${repo_root}/models/manifest.yaml" "${bundle}/model-manifest.yaml"
cp "${repo_root}/references/manifest.yaml" "${bundle}/reference-manifest.yaml"
cp -R "${repo_root}/data/manifests" "${bundle}/data-manifests"
cp "${repo_root}/src/cas13_if/refold/output_schema.json" "${bundle}/schemas/"

if [[ -n "${input_root}" ]]; then
  if [[ ! -d "${input_root}" ]]; then
    echo "ERROR: input shard directory does not exist: ${input_root}" >&2
    exit 4
  fi
  while IFS= read -r -d '' input_path; do
    relative="${input_path#"${input_root}/"}"
    mkdir -p "${bundle}/inputs/$(dirname "${relative}")"
    cp "${input_path}" "${bundle}/inputs/${relative}"
  done < <(
    find "${input_root}" -type f \
      \( -name '*.fasta' -o -name '*.fa' -o -name '*.json' -o -name '*.jsonl' \
      -o -name '*.yaml' -o -name '*.yml' \) -print0
  )
fi

printf '%s\n' \
  "git clone https://github.com/XuJP264/ESM-IF-on-cas13.git ESM-IF" \
  "cd ESM-IF" \
  "git checkout ${commit}" \
  > "${bundle}/repo-clone.txt"

python3 - "${repo_root}" "${bundle}" "${commit}" "${config_hash}" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
bundle = Path(sys.argv[2])
commit = sys.argv[3]
config_hash = sys.argv[4]
asset_paths = (
    "models/esm_if1/esm_if1_gvp4_t16_142M_UR50.pt",
    "models/proteinmpnn/v_48_020.pt",
    "models/ligandmpnn/proteinmpnn_v_48_020.pt",
    "models/ligandmpnn/ligandmpnn_v_32_010_25.pt",
    "models/ligandmpnn/solublempnn_v_48_020.pt",
    "data/raw/atlas/v1.0/crispr-cas-atlas-v1.0.json",
    "data/experimental_structures/6e9f.pdb",
    "data/experimental_structures/6e9f.cif",
    "data/experimental_structures/5xwp.pdb",
    "data/experimental_structures/5xwp.cif",
    "data/experimental_structures/6e9e.pdb",
    "data/experimental_structures/6e9e.cif",
    "data/experimental_structures/5xwy.pdb",
    "data/experimental_structures/5xwy.cif",
)
assets = []
missing = []
for relative in asset_paths:
    path = root / relative
    if path.is_file():
        import hashlib

        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        assets.append(
            {
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": digest.hexdigest(),
            }
        )
    else:
        missing.append(relative)
dirty = bool(
    __import__("subprocess")
    .run(
        ["git", "-C", str(root), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    )
    .stdout.strip()
)
manifest = {
    "schema_version": "1.0",
    "git_commit": commit,
    "git_worktree_dirty_at_export": dirty,
    "config_hash": config_hash,
    "is_mock": False,
    "large_assets_embedded": False,
    "assets": assets,
    "missing_assets": missing,
    "input_files": sorted(
        str(path.relative_to(bundle))
        for path in (bundle / "inputs").rglob("*")
        if path.is_file()
    ),
    "transfer": (
        "Transfer this small bundle separately; use scripts/sync_assets.sh "
        "for large assets and verify their SHA256 values."
    ),
}
(bundle / "bundle-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
(bundle / "ASSET_SHA256SUMS").write_text(
    "".join(f"{asset['sha256']}  {asset['path']}\n" for asset in assets)
)
PY
(cd "${bundle}" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
echo "${bundle}"
