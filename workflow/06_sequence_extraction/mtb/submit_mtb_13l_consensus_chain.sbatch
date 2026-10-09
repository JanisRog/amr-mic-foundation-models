#!/usr/bin/env bash
#SBATCH --job-name=mtb_13l_chain
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
PIPE="${PIPE:?set PIPE in config/paths.env}"
PYTHON_BIN="${PYTHON_BIN:-${PY:?set PY in config/paths.env}}"
SNIPPY_SIF="${SNIPPY_SIF:?set SNIPPY_SIF in config/paths.env}"
REF="${REF:-$MTB/data/reference/tbprofiler_tbdb_72ef6fa/tbdb.fasta}"

SNIPPY_RUNS="${SNIPPY_RUNS:-$MTB/results/gwas/snippy_h37rv/runs}"
CLUMPS="${CLUMPS:-$MTB/results/gwas/selections/mtb_all_gene_target_clumps/all_gene_target_clumps.tsv}"
SAMPLE_LIST="${SAMPLE_LIST:-$MTB/results/gwas/snippy_h37rv/qc/ok_samples.txt}"
OUTROOT="${OUTROOT:-$MTB/results/gwas/sample_windows/mtb_all_gene_consensus_chain}"

CONSENSUS_NAME="${CONSENSUS_NAME:-snps.consensus.fa}"
CHAIN_NAME="${CHAIN_NAME:-snps.consensus.chain}"
VCF_GZ_NAME="${VCF_GZ_NAME:-snps.vcf.gz}"
VCF_NAME="${VCF_NAME:-snps.vcf}"
REGENERATE="${REGENERATE:-1}"

SCRIPT="$PIPE/scripts/13l_build_snp_gene_consensus_chain_sequences.py"

for f in "$SCRIPT" "$REF" "$CLUMPS" "$SAMPLE_LIST" "$SNIPPY_SIF"; do
  [[ -s "$f" ]] || { echo "[ERROR] missing: $f" >&2; exit 2; }
done
[[ -d "$SNIPPY_RUNS" ]] || { echo "[ERROR] missing snippy runs dir: $SNIPPY_RUNS" >&2; exit 2; }

mkdir -p "$OUTROOT" "$MTB/00_logs/mtb_13l"
STAMP="$(date -u +%Y%m%d_%H%M%S)"
RUN_OUT="$OUTROOT/run_${STAMP}"
mkdir -p "$RUN_OUT"

echo "[INFO] clumps        = $CLUMPS ($(( $(wc -l < "$CLUMPS") - 1 )) targets)"
echo "[INFO] sample_list   = $SAMPLE_LIST ($(grep -cve '^[[:space:]]*$' "$SAMPLE_LIST") samples)"
echo "[INFO] sample_root   = $SNIPPY_RUNS"
echo "[INFO] reference     = $REF"
echo "[INFO] run_out       = $RUN_OUT"

# 13l shells out to bcftools; route it through the snippy container.
module load apptainer >/dev/null 2>&1 || true
CTR_BIN="$(command -v apptainer || command -v singularity || true)"
[[ -n "$CTR_BIN" ]] || { echo "[ERROR] no apptainer/singularity on PATH after 'module load apptainer'" >&2; exit 2; }
echo "[INFO] container_bin = $CTR_BIN"
BCFTOOLS_WRAPPER="$RUN_OUT/bcftools_wrapper.sh"
cat > "$BCFTOOLS_WRAPPER" <<BCFEOF
#!/usr/bin/env bash
exec "$CTR_BIN" exec --cleanenv -B "${BIND_PATHS:?set BIND_PATHS in config/paths.env}" "$SNIPPY_SIF" bcftools "\$@"
BCFEOF
chmod +x "$BCFTOOLS_WRAPPER"

REGEN_FLAG=()
[[ "$REGENERATE" == "1" ]] && REGEN_FLAG=(--regenerate-consensus)

"$PYTHON_BIN" "$SCRIPT" \
  --target-clumps "$CLUMPS" \
  --sample-root "$SNIPPY_RUNS" \
  --reference-fasta "$REF" \
  --sample-ids "$SAMPLE_LIST" \
  --consensus-name "$CONSENSUS_NAME" \
  --chain-name "$CHAIN_NAME" \
  --sample-vcf-gz-name "$VCF_GZ_NAME" \
  --sample-vcf-name "$VCF_NAME" \
  --bcftools-bin "$BCFTOOLS_WRAPPER" \
  --outdir "$RUN_OUT" \
  "${REGEN_FLAG[@]}"

ln -sfn "$RUN_OUT" "$OUTROOT/latest"
echo "[OK] mtb_13l_latest=$OUTROOT/latest"
