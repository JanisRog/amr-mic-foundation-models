#!/usr/bin/env bash
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GWAS_PIPELINE_ROOT="$(cd "$THIS_DIR/.." && pwd)"
AMR_ROOT="$(cd "$GWAS_PIPELINE_ROOT/.." && pwd)"

CONFIG_ENV="${CONFIG_ENV:-}"
[[ -n "$CONFIG_ENV" ]] || { echo "[ERROR] CONFIG_ENV not set" >&2; exit 1; }
[[ -f "$CONFIG_ENV" ]] || { echo "[ERROR] CONFIG_ENV missing: $CONFIG_ENV" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG_ENV"

require_file() { [[ -f "$1" ]] || { echo "[ERROR] missing file: $1" >&2; exit 1; }; }
require_dir() { [[ -d "$1" ]] || { echo "[ERROR] missing directory: $1" >&2; exit 1; }; }

refresh_symlink() {
  local target="$1"
  local link_path="$2"

  mkdir -p "$(dirname "$link_path")"

  if [[ -L "$link_path" || -f "$link_path" ]]; then
    ln -sfn "$target" "$link_path"
    return 0
  fi

  if [[ -d "$link_path" ]]; then
    local backup="${link_path}.bak_$(date +%Y%m%d_%H%M%S)"
    mv "$link_path" "$backup"
    echo "[WARN] moved existing directory aside to refresh symlink: $link_path -> $backup"
  fi

  ln -sfn "$target" "$link_path"
}

resolve_container_bin() {
  if command -v apptainer >/dev/null 2>&1; then
    echo "$(command -v apptainer)"
    return 0
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo "$(command -v singularity)"
    return 0
  fi
  module load apptainer >/dev/null 2>&1 || module load Apptainer >/dev/null 2>&1 || module load singularityce >/dev/null 2>&1 || true
  if command -v apptainer >/dev/null 2>&1; then
    echo "$(command -v apptainer)"
    return 0
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo "$(command -v singularity)"
    return 0
  fi
  echo "[ERROR] neither apptainer nor singularity found" >&2
  exit 1
}

setup_stage_logging() {
  local stage="$1"
  local log_dir="${DATASET_LOG_ROOT}/gwas"
  mkdir -p "$log_dir" "$GWAS_PIPELINE_ROOT/logs" "$RESULTS_ROOT" "$GWAS_ROOT" "$MODELS_ROOT" "$DEPLOY_ROOT"
  local stem="${stage}_${SLURM_JOB_ID:-$$}"
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    stem+="_${SLURM_ARRAY_TASK_ID}"
  fi
  exec >>"$log_dir/${stem}.out" 2>>"$log_dir/${stem}.err"
  echo "[INFO] stage=$stage"
  echo "[INFO] timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "[INFO] dataset_name=$DATASET_NAME"
  echo "[INFO] dataset_root=$DATASET_ROOT"
  echo "[INFO] config_env=$CONFIG_ENV"
}
