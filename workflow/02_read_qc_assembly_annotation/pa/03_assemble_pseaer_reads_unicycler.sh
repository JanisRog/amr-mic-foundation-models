#!/usr/bin/env bash
#SBATCH --job-name=pseaer_assembly
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G
#SBATCH --time=24:00:00
#SBATCH --output=03_assemble_pseaer_%A_%a.out
#SBATCH --error=03_assemble_pseaer_%A_%a.err

set -euo pipefail

BASE="${BASE:-${PSEAER:?set PSEAER in config/paths.env}}"
REPRO="${REPRO:-$BASE/reproducibility}"
ACC_FILE="${ACC_FILE:-$BASE/metadata/parsed/accessions_all_direct.txt}"
RAW_ROOT="${RAW_ROOT:-$BASE/raw_reads}"
FASTQ_ROOT="${FASTQ_ROOT:-$RAW_ROOT/fastq}"
OUT_ROOT="${OUT_ROOT:-$BASE/assemblies/unicycler_pilon}"
STATUS_DIR="${STATUS_DIR:-$BASE/assemblies/status}"
TMP_ROOT="${TMP_ROOT:-$BASE/assemblies/tmp}"
FORCE="${FORCE:-0}"

SING_ROOT="${SING_ROOT:?set SING_ROOT in config/paths.env}"
TRIMMOMATIC_IMG="${TRIMMOMATIC_IMG:-$SING_ROOT/trimmomatic-0.39--hdfd78af_2}"
UNICYCLER_IMG="${UNICYCLER_IMG:-$SING_ROOT/unicycler_0.5.0}"
BWA_IMG="${BWA_IMG:-$SING_ROOT/quay.io-biocontainers-bwa-0.7.17--h5bf99c6_8}"
SAMTOOLS_IMG="${SAMTOOLS_IMG:-$SING_ROOT/samtools_1.15.1}"
PILON_IMG="${PILON_IMG:-$SING_ROOT/quay.io-biocontainers-pilon-1.24--hdfd78af_0}"

ILLUMINACLIP="${ILLUMINACLIP:?set ILLUMINACLIP in config/paths.env}"
TRIMMOMATIC_PE_EXTRA="${TRIMMOMATIC_PE_EXTRA:-SLIDINGWINDOW:4:12 MINLEN:100}"
TRIMMOMATIC_SE_EXTRA="${TRIMMOMATIC_SE_EXTRA:-SLIDINGWINDOW:4:12 MINLEN:100}"
PILON_JAVA_XMX="${PILON_JAVA_XMX:-10g}"

mkdir -p "$OUT_ROOT" "$STATUS_DIR" "$TMP_ROOT" "$REPRO/logs"

[[ -f "$ACC_FILE" ]] || { echo "[ERROR] missing accession file: $ACC_FILE" >&2; exit 1; }
if [[ -n "${ACC_ID:-}" ]]; then
  ACC="$ACC_ID"
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  ACC="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACC_FILE" | tr -d '\r')"
else
  echo "[ERROR] set ACC_ID or submit as a Slurm array with ACC_FILE" >&2
  exit 1
fi
[[ -n "$ACC" ]] || { echo "[ERROR] empty accession at task id ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

OK_MARK="$STATUS_DIR/${ACC}.ok"
FAIL_MARK="$STATUS_DIR/${ACC}.fail"
SAMPLE_FASTQ_DIR="$FASTQ_ROOT/$ACC"
SAMPLE_OUT="$OUT_ROOT/$ACC"
SAMPLE_TMP="$TMP_ROOT/${ACC}_${SLURM_JOB_ID:-$$}"
WORK_DIR="$SAMPLE_TMP/work"
ASSEMBLY_FINAL="$SAMPLE_OUT/${ACC}_assembly.fasta"

RAW_R1="$SAMPLE_FASTQ_DIR/${ACC}_1.fastq.gz"
RAW_R2="$SAMPLE_FASTQ_DIR/${ACC}_2.fastq.gz"
RAW_SE="$SAMPLE_FASTQ_DIR/${ACC}.fastq.gz"

cleanup_fail() {
  touch "$FAIL_MARK"
}
trap cleanup_fail ERR

if [[ "$FORCE" != "1" && -s "$ASSEMBLY_FINAL" ]]; then
  echo "[OK] accession=$ACC already assembled"
  touch "$OK_MARK"
  exit 0
fi

module purge || true
if ! command -v apptainer >/dev/null 2>&1; then
  module load apptainer >/dev/null 2>&1 || module load Apptainer >/dev/null 2>&1 || true
fi
command -v apptainer >/dev/null 2>&1 || { echo "[ERROR] apptainer not found in PATH" >&2; exit 1; }

unset SINGULARITY_BINDPATH || true
export APPTAINER_BINDPATH="${BIND_PATHS:?set BIND_PATHS in config/paths.env}"
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
export APPTAINER_CACHEDIR="${HOME}/.apptainer/cache"
mkdir -p "$APPTAINER_CACHEDIR"

for img in "$TRIMMOMATIC_IMG" "$UNICYCLER_IMG" "$BWA_IMG" "$SAMTOOLS_IMG" "$PILON_IMG"; do
  [[ -e "$img" ]] || { echo "[ERROR] missing container: $img" >&2; exit 1; }
done
[[ -f "$ILLUMINACLIP" ]] || { echo "[ERROR] missing adapter file: $ILLUMINACLIP" >&2; exit 1; }

if [[ "$FORCE" != "1" && -d "$SAMPLE_OUT" && ! -s "$ASSEMBLY_FINAL" ]]; then
  echo "[INFO] accession=$ACC has existing incomplete output dir; proceeding with fresh temp work dir"
fi

rm -rf "$WORK_DIR"
mkdir -p "$SAMPLE_OUT" "$SAMPLE_TMP" "$WORK_DIR"
cd "$WORK_DIR"
rm -f "$FAIL_MARK"

echo "[INFO] accession=$ACC"
echo "[INFO] sample_fastq_dir=$SAMPLE_FASTQ_DIR"
echo "[INFO] sample_out=$SAMPLE_OUT"
echo "[INFO] sample_tmp=$SAMPLE_TMP"

run_in_img() {
  local img="$1"
  shift
  apptainer exec --bind "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$img" "$@"
}

if [[ -f "$RAW_R1" && -f "$RAW_R2" ]]; then
  mode="paired"
elif [[ -f "$RAW_SE" ]]; then
  mode="single"
else
  echo "[ERROR] no expected FASTQ files found for $ACC in $SAMPLE_FASTQ_DIR" >&2
  exit 1
fi

if [[ "$mode" == "paired" ]]; then
  trim1p="${ACC}_R1_trim_paired.fastq.gz"
  trim2p="${ACC}_R2_trim_paired.fastq.gz"
  trim1u="${ACC}_R1_trim_unpaired.fastq.gz"
  trim2u="${ACC}_R2_trim_unpaired.fastq.gz"

  echo "[INFO] trimming paired reads $ACC $(date +\"%T\")"
  run_in_img "$TRIMMOMATIC_IMG" trimmomatic PE \
    -threads "$SLURM_CPUS_PER_TASK" \
    -phred33 \
    "$RAW_R1" "$RAW_R2" \
    "$trim1p" "$trim1u" "$trim2p" "$trim2u" \
    "ILLUMINACLIP:${ILLUMINACLIP}:2:30:10:2" \
    $TRIMMOMATIC_PE_EXTRA \
    2> "${ACC}.quality_read_trimm_info"

  echo "[INFO] assembling paired reads $ACC $(date +\"%T\")"
  run_in_img "$UNICYCLER_IMG" unicycler \
    -1 "$trim1p" \
    -2 "$trim2p" \
    -t "$SLURM_CPUS_PER_TASK" \
    -o ./

  ASSEMBLY_IN="assembly.fasta"
  ASSEMBLY_DRAFT="${ACC}_draft_assembly.fasta"
  mv "$ASSEMBLY_IN" "$ASSEMBLY_DRAFT"
  [[ -f unicycler.log ]] && mv unicycler.log "${ACC}_unicycler.log"

  echo "[INFO] remapping paired reads $ACC $(date +\"%T\")"
  run_in_img "$BWA_IMG" bwa index "$ASSEMBLY_DRAFT"
  run_in_img "$BWA_IMG" bwa mem -t "$SLURM_CPUS_PER_TASK" "$ASSEMBLY_DRAFT" "$trim1p" "$trim2p" > "${ACC}_alignment.sam"
else
  trim1="${ACC}_trimmed.fastq.gz"

  echo "[INFO] trimming single-end reads $ACC $(date +\"%T\")"
  run_in_img "$TRIMMOMATIC_IMG" trimmomatic SE \
    -threads "$SLURM_CPUS_PER_TASK" \
    -phred33 \
    "$RAW_SE" "$trim1" \
    "ILLUMINACLIP:${ILLUMINACLIP}:2:30:10" \
    $TRIMMOMATIC_SE_EXTRA \
    2> "${ACC}.quality_read_trimm_info"

  echo "[INFO] assembling single-end reads $ACC $(date +\"%T\")"
  run_in_img "$UNICYCLER_IMG" unicycler \
    -s "$trim1" \
    -t "$SLURM_CPUS_PER_TASK" \
    -o ./

  ASSEMBLY_IN="assembly.fasta"
  ASSEMBLY_DRAFT="${ACC}_draft_assembly.fasta"
  mv "$ASSEMBLY_IN" "$ASSEMBLY_DRAFT"
  [[ -f unicycler.log ]] && mv unicycler.log "${ACC}_unicycler.log"

  echo "[INFO] remapping single-end reads $ACC $(date +\"%T\")"
  run_in_img "$BWA_IMG" bwa index "$ASSEMBLY_DRAFT"
  run_in_img "$BWA_IMG" bwa mem -t "$SLURM_CPUS_PER_TASK" "$ASSEMBLY_DRAFT" "$trim1" > "${ACC}_alignment.sam"
fi

run_in_img "$SAMTOOLS_IMG" samtools sort -@ "$SLURM_CPUS_PER_TASK" -T "${ACC}" -o "${ACC}_alignment.bam" "${ACC}_alignment.sam"
run_in_img "$SAMTOOLS_IMG" samtools rmdup "${ACC}_alignment.bam" "${ACC}_alignment.removed_duplicates.bam"
run_in_img "$SAMTOOLS_IMG" samtools index "${ACC}_alignment.removed_duplicates.bam"

echo "[INFO] pilon polishing $ACC $(date +\"%T\")"
run_in_img "$PILON_IMG" /bin/bash -lc "export _JAVA_OPTIONS='-Xmx${PILON_JAVA_XMX}'; pilon --mindepth 5 --minmq 10 --threads ${SLURM_CPUS_PER_TASK} --genome '${ASSEMBLY_DRAFT}' --frags '${ACC}_alignment.removed_duplicates.bam' --changes --variant --outdir . --output '${ACC}_assembly'"

[[ -f "${ACC}_assembly.fasta" ]] || { echo "[ERROR] Pilon assembly missing for $ACC" >&2; exit 1; }

cp -a "$WORK_DIR"/. "$SAMPLE_OUT"/
[[ -f "$ASSEMBLY_FINAL" ]] || { echo "[ERROR] final assembly missing after sync for $ACC" >&2; exit 1; }

rm -rf "$SAMPLE_TMP"
rm -f "$FAIL_MARK"
touch "$OK_MARK"

echo "[OK] accession=$ACC done $(date +\"%T\")"
echo "[OK] assembly=$ASSEMBLY_FINAL"
