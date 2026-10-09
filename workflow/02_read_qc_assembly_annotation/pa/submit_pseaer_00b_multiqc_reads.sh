#!/usr/bin/env bash
#SBATCH --job-name=pseaer_multiqc
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail
CONFIG_ENV="${PIPE:?set PIPE in config/paths.env}/00_submit/pseaer.env"
WRAPPER="${PIPE:?set PIPE in config/paths.env}/wrappers/00b_multiqc.sh"
export CONFIG_ENV
exec bash "$WRAPPER"
