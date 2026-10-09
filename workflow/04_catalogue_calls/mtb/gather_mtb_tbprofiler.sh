#!/usr/bin/env bash
#SBATCH --job-name=mtb_tbprof_gather
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

MTB_ROOT="${MTB_ROOT:-${MTB:?set MTB in config/paths.env}}"
OUTROOT="${OUTROOT:-$MTB_ROOT/results/tbprofiler/reads}"
TBPROFILER_SIF="${TBPROFILER_SIF:?set TBPROFILER_SIF in config/paths.env}"
SAMPLE_LIST="${SAMPLE_LIST:-$MTB_ROOT/reads/paired_sample_ids.txt}"
CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-}"

resolve_container_runtime() {
  if command -v apptainer >/dev/null 2>&1; then command -v apptainer; return 0; fi
  if command -v singularity >/dev/null 2>&1; then command -v singularity; return 0; fi
  module load apptainer >/dev/null 2>&1 || module load Apptainer >/dev/null 2>&1 || module load singularityce >/dev/null 2>&1 || true
  if command -v apptainer >/dev/null 2>&1; then command -v apptainer; return 0; fi
  if command -v singularity >/dev/null 2>&1; then command -v singularity; return 0; fi
  return 1
}

if [[ -z "$CONTAINER_RUNTIME" ]]; then
  if ! CONTAINER_RUNTIME="$(resolve_container_runtime)"; then
    echo "[ERROR] neither apptainer nor singularity found; set CONTAINER_RUNTIME to its absolute path" >&2
    exit 2
  fi
fi

mkdir -p "$OUTROOT/gather/results"
find "$OUTROOT/samples" -path '*/results/*.results.json' -type f -print0 |
while IFS= read -r -d '' path; do
  ln -sfn "$path" "$OUTROOT/gather/results/$(basename "$path")"
done

expected=$(wc -l < "$SAMPLE_LIST")
complete=$(find "$OUTROOT/gather/results" -name '*.results.json' -type l | wc -l)
echo "[INFO] expected=$expected complete=$complete missing=$((expected - complete))"

cd "$OUTROOT/gather"
unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
"$CONTAINER_RUNTIME" exec --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$TBPROFILER_SIF" tb-profiler collate --dir results --prefix tbprofiler
echo "[OK] gathered=$OUTROOT/gather"
