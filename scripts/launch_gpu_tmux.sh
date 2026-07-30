#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 CONFIG TASK" >&2
  exit 2
fi
if ! command -v tmux >/dev/null 2>&1; then
  echo "ERROR: tmux is required" >&2
  exit 3
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
config="$1"
task="$2"
timestamp="$(date +%Y%m%d-%H%M%S)"
safe_task="$(printf '%s' "${task}" | tr -c 'A-Za-z0-9_-' '-')"
session="cas13-${safe_task}-${timestamp}"
run_dir="${repo_root}/results/runs/${timestamp}-${safe_task}-gpu"
mkdir -p "${run_dir}"
printf '%s\n' "${session}" > "${repo_root}/results/latest_gpu_session.txt"
printf '%s\n' "${run_dir}" > "${repo_root}/results/latest_gpu_run.txt"

command="cd $(printf '%q' "${repo_root}") && set -o pipefail && make $(printf '%q' "${task}") CONFIG=$(printf '%q' "${config}") > $(printf '%q' "${run_dir}/stdout.log") 2> $(printf '%q' "${run_dir}/stderr.log"); code=\$?; printf '%s\n' \"\$code\" > $(printf '%q' "${run_dir}/EXIT_CODE"); exit \"\$code\""
tmux new-session -d -s "${session}" "${command}"
echo "session=${session}"
echo "run_dir=${run_dir}"

