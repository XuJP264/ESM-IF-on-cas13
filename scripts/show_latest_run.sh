#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
find "${repo_root}/results/runs" -mindepth 1 -maxdepth 1 -type d \
  -printf '%T@ %p\n' | sort -nr | head -n 1 | cut -d' ' -f2-

