#!/usr/bin/env bash
#SBATCH --job-name=pseaer_bakta_ref
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00

set -euo pipefail
CONFIG_ENV="${PIPE:?set PIPE in config/paths.env}/00_submit/pseaer.env"
WRAPPER="${PIPE:?set PIPE in config/paths.env}/wrappers/02r_bakta_reference.sh"
export CONFIG_ENV
exec bash "$WRAPPER"
