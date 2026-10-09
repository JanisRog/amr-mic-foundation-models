#!/usr/bin/env bash
#SBATCH --job-name=pseaer_quast
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00

set -euo pipefail
CONFIG_ENV="${PIPE:?set PIPE in config/paths.env}/00_submit/pseaer.env"
WRAPPER="${PIPE:?set PIPE in config/paths.env}/wrappers/00c_quast.sh"
export CONFIG_ENV
exec bash "$WRAPPER"
