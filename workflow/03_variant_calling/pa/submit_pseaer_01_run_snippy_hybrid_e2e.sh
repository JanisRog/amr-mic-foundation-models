#!/usr/bin/env bash
#SBATCH --job-name=pseaer_snippy
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00

set -euo pipefail
CONFIG_ENV="${PIPE:?set PIPE in config/paths.env}/00_submit/pseaer.env"
WRAPPER="${PIPE:?set PIPE in config/paths.env}/wrappers/01_run_snippy_hybrid_e2e.sh"
export CONFIG_ENV
exec bash "$WRAPPER"
