#!/usr/bin/env bash
# Build a cluster-side analysis root for the MTB holdout sweeps.
#
# Run this ON the cluster. The big inputs (embeddings, seqid maps) stay where
# they were produced and are symlinked in; only the small files -- splits,
# configs, metadata, tbprofiler binary matrices, and the intergenic-only GWAS
# map -- need uploading from the Mac.
set -euo pipefail

MTB="${MTB:?set MTB in config/paths.env}"
ROOT="${ROOT:-$MTB/results/ml/03_mtb}"
EMB_RUN="${EMB_RUN:-$MTB/results/gwas/carbon_embeddings/carbon_8b/mtb_all_gene_consensus_chain/latest}"
CHAIN="${CHAIN:-$MTB/results/gwas/sample_windows/mtb_all_gene_consensus_chain/latest}"
GWAS_EMB_RUN="${GWAS_EMB_RUN-$MTB/results/gwas/carbon_embeddings/carbon_8b/plos_mic_remaining_gene_up300_intergenic/run_20260810_130033}"
GWAS_CHAIN="${GWAS_CHAIN:-$MTB/results/gwas/sample_windows/plos_mic_remaining_gene_up300_intergenic/latest}"
ANTIBIOTIC="${ANTIBIOTIC:-rifampicin}"
PYTHON_BIN="${PYTHON_BIN:-${PY:?set PY in config/paths.env}}"
# rows of each map to check; 0 = all (slow on a login node)
VALIDATE_ROWS="${VALIDATE_ROWS:-2000000}"
# Borrow the small model-independent inputs (metadata, binary matrices, GWAS map,
# splits) from an existing root instead of re-uploading them per model.
SHARE_FROM="${SHARE_FROM:-}"
ANTIBIOTICS="${ANTIBIOTICS:-rifampicin isoniazid ethambutol}"

# An embedding root may be given as the family directory rather than the run
# inside it; accept either by descending into latest/ when needed.
resolve_run() {
  local d="$1"
  [[ -n "$d" ]] || { printf ''; return 0; }
  if [[ ! -s "$d/embeddings.npy" ]]; then
    if [[ -s "$d/latest/embeddings.npy" ]]; then
      d="$d/latest"
    else
      # not every family dir has a latest/ symlink; fall back to the newest run_*
      local newest
      newest="$(find -L "$d" -maxdepth 1 -type d -name 'run_*' 2>/dev/null | sort | tail -1)"
      [[ -n "$newest" && -s "$newest/embeddings.npy" ]] && d="$newest"
    fi
  fi
  printf '%s' "$d"
}
EMB_RUN="$(resolve_run "$EMB_RUN")"
GWAS_EMB_RUN="$(resolve_run "$GWAS_EMB_RUN")"

EMB_DST="$ROOT/ml_inputs/embeddings/tb_profiler"
GWAS_DST="$ROOT/ml_inputs/embeddings/gwas"
mkdir -p "$EMB_DST" "$GWAS_DST" \
         "$ROOT/ml_inputs/tbprofiler_binary" "$ROOT/metadata" \
         "$ROOT/$ANTIBIOTIC/results/holdout_results"

link() {  # link <src> <dst>
  [[ -s "$1" ]] || { echo "[ERROR] missing source: $1" >&2; exit 2; }
  ln -sfn "$1" "$2"
  echo "[link] $(basename "$2") -> $1"
}

link "$EMB_RUN/embeddings.npy"            "$EMB_DST/embeddings.npy"
link "$EMB_RUN/embedding_index.tsv"       "$EMB_DST/embedding_index.tsv"
link "$CHAIN/sample_target_seqid_map.tsv" "$EMB_DST/sample_target_seqid_map.tsv"

# evo2 halves sequences over 8192 bp; the sweep reassembles them via this map,
# so it must sit beside the embeddings or those targets cannot resolve.
# It is written by the splitter, so it may live either with the embeddings or
# with the chain windows it was derived from -- look in both, resolving symlinks
# so a "latest" pointer does not hide it.
# Explicit map paths, for a model that embedded someone else's split fasta and so
# has no map of its own beside its embeddings (DNABERT-2 reuses evo2's halves).
SPLIT_MAP="${SPLIT_MAP:-}"
GWAS_SPLIT_MAP="${GWAS_SPLIT_MAP:-}"

link_split_map() {  # link_split_map <dst_dir> <search_dir>...
  local dst="$1"; shift
  local dirs=() d sm
  local EXPLICIT_MAP="${EXPLICIT_MAP:-}"
  for d in "$@"; do [[ -n "$d" && -d "$d" ]] && dirs+=("$d"); done
  sm=""
  if [[ -n "$EXPLICIT_MAP" ]]; then
    [[ -s "$EXPLICIT_MAP" ]] || { echo "[ERROR] SPLIT_MAP not readable: $EXPLICIT_MAP" >&2; exit 2; }
    sm="$EXPLICIT_MAP"
    dirs=()
  fi
  if [[ ${#dirs[@]} -gt 0 ]]; then
    sm="$(find -L "${dirs[@]}" -maxdepth 1 \
          -name 'unique_window_sequences.split_long_*.map.tsv' 2>/dev/null | head -1)"
  fi
  # a stale map from a previous, differently-split model would silently
  # mis-resolve, so clear before linking
  rm -f "$dst"/unique_window_sequences.split_long_*.map.tsv
  [[ -n "$sm" ]] && link "$sm" "$dst/$(basename "$sm")"
  return 0
}

EXPLICIT_MAP="$SPLIT_MAP"      link_split_map "$EMB_DST"  "$EMB_RUN" "$CHAIN"

# GWAS block: Carbon embeddings for the PLOS gene/upstream/intergenic windows.
# The map is uploaded separately (intergenic-only) so it is never overwritten.
if [[ -n "$GWAS_EMB_RUN" && -s "$GWAS_EMB_RUN/embeddings.npy" ]]; then
  link "$GWAS_EMB_RUN/embeddings.npy"      "$GWAS_DST/embeddings.npy"
  link "$GWAS_EMB_RUN/embedding_index.tsv" "$GWAS_DST/embedding_index.tsv"
  # this block is split for evo2 too, so it needs its own halves map
  EXPLICIT_MAP="$GWAS_SPLIT_MAP" link_split_map "$GWAS_DST" "$GWAS_EMB_RUN" "$GWAS_CHAIN"
else
  EXPLICIT_MAP="" link_split_map "$GWAS_DST"
  # Clear any links a previous build left here. Otherwise disabling the block
  # silently keeps the *previous model's* embeddings wired in, and the dims
  # match (4096 for both Carbon and evo2) so nothing downstream complains.
  for f in embeddings.npy embedding_index.tsv; do
    if [[ -e "$GWAS_DST/$f" || -L "$GWAS_DST/$f" ]]; then
      rm -f "$GWAS_DST/$f"; echo "[clear] gwas/$f"
    fi
  done
  echo "[skip] gwas embeddings not wired (GWAS_EMB_RUN=$GWAS_EMB_RUN)"
fi

# --- share small inputs from an existing root ------------------------------
if [[ -n "$SHARE_FROM" ]]; then
  [[ -d "$SHARE_FROM" ]] || { echo "[ERROR] SHARE_FROM not a directory: $SHARE_FROM" >&2; exit 2; }
  for rel in metadata ml_inputs/tbprofiler_binary; do
    if [[ -e "$SHARE_FROM/$rel" ]]; then
      [[ -e "$ROOT/$rel" || -L "$ROOT/$rel" ]] && rm -rf "$ROOT/$rel"
      mkdir -p "$(dirname "$ROOT/$rel")"
      ln -sfn "$SHARE_FROM/$rel" "$ROOT/$rel"
      echo "[share] $rel"
    fi
  done
  gmap="$SHARE_FROM/ml_inputs/embeddings/gwas/sample_target_seqid_map.tsv"
  [[ -e "$gmap" ]] && ln -sfn "$(readlink -f "$gmap")" "$GWAS_DST/sample_target_seqid_map.tsv" \
    && echo "[share] gwas/sample_target_seqid_map.tsv"
  for ab in $ANTIBIOTICS; do
    src="$SHARE_FROM/$ab/results/holdout_results"
    dst="$ROOT/$ab/results/holdout_results"
    [[ -d "$src" ]] || continue
    mkdir -p "$dst"
    n=0
    while IFS= read -r f; do
      ln -sfn "$(readlink -f "$f")" "$dst/$(basename "$f")"; n=$((n+1))
    done < <(find "$src" -maxdepth 1 \( -type f -o -type l \) \
             \( -name 'sample_splits_*' -o -name 'sampled_configs_*' \))
    echo "[share] $ab: $n split/config files"
  done
fi

[[ -x "$PYTHON_BIN" ]] || { echo "[ERROR] PYTHON_BIN not executable: $PYTHON_BIN" >&2; exit 2; }

cat > "$ROOT/.validate_block.py" <<'PYEOF'
import sys, csv, pathlib, numpy as np

d = pathlib.Path(sys.argv[1])
limit = int(sys.argv[2])


def npy_shape(path):
    """Read shape/dtype from the .npy header only.

    mmap of a multi-GB array fails on the login node's address-space cap, and
    the header is all this check needs.
    """
    with open(path, "rb") as fh:
        ver = np.lib.format.read_magic(fh)
        if ver == (1, 0):
            shape, _, dtype = np.lib.format.read_array_header_1_0(fh)
        elif ver == (2, 0):
            shape, _, dtype = np.lib.format.read_array_header_2_0(fh)
        else:
            raise RuntimeError("unsupported .npy version %s" % (ver,))
    return shape, dtype


emb_shape, emb_dtype = npy_shape(d / "embeddings.npy")
ids = set()
with (d / "embedding_index.tsv").open() as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        ids.add(row["seq_id"])

# sequences split into parts resolve through the split map, not directly.
# Column names have varied between splitter versions, so pick them by content
# rather than by a hardcoded pair, and say which pair was used.
ORIG_NAMES = ("original_sequence_id", "orig_seq_id", "source_seq_id", "seq_id")
PART_NAMES = ("split_sequence_id", "part_seq_id", "chunk_seq_id", "split_id")


def pick(header, names, exclude=()):
    for n in names:
        if n in header and n not in exclude:
            return n
    return None


via_split = set()
for sm in d.glob("unique_window_sequences.split_long_*.map.tsv"):
    with sm.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        header = rdr.fieldnames or []
        part_col = pick(header, PART_NAMES)
        orig_col = pick(header, ORIG_NAMES, exclude=(part_col,))
        if not part_col or not orig_col:
            print(f"[ERROR] {sm.name}: cannot identify columns in {header}")
            sys.exit(1)
        parts = {}
        for row in rdr:
            parts.setdefault(row[orig_col], []).append(row[part_col])
    for orig, ps in parts.items():
        if orig and all(x in ids for x in ps):
            via_split.add(orig)
    print(f"[OK] split map {sm.name}: {len(parts):,} split sequences "
          f"({orig_col} -> {part_col}), {len(via_split):,} fully embedded")
if not via_split:
    print("[WARN] no split map linked -- long sequences will not resolve")
resolvable = ids | via_split
print(f"[OK] embeddings shape={emb_shape} dtype={emb_dtype}")
print(f"[OK] embedding_index seq_ids={len(ids):,}")
if len(ids) != emb_shape[0]:
    print(f"[ERROR] index rows ({len(ids)}) != embedding rows ({emb_shape[0]})")
    sys.exit(1)

# stream the map; it may be far too large to hold in memory here
seen, missing, rows, targets = set(), 0, 0, set()
bad = {}
with (d / "sample_target_seqid_map.tsv").open() as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        rows += 1
        targets.add(row.get("target_id", ""))
        s = row.get("seq_id", "")
        if s in resolvable:
            seen.add(s)
        else:
            missing += 1
            if len(bad) < 5:
                bad.setdefault(s, row.get("target_id", ""))
        if limit and rows >= limit:
            break
scope = "all" if not limit else f"first {limit:,}"
print(f"[OK] map rows checked={rows:,} ({scope})  targets={len(targets):,}  "
      f"distinct seq_ids={len(seen):,}  unresolved={missing:,}")
if missing:
    print(f"[ERROR] map references sequences absent from the embedding index "
          f"({len(set(bad)):,}+ distinct); examples:")
    for s, t in bad.items():
        print(f"        seq_id={s}  target_id={t}")
    sys.exit(1)
print("[OK] every mapped sequence resolves to an embedding row")
PYEOF

validate_block() {  # validate_block <dir> <label>
  local dir="$1" label="$2"
  if [[ ! -s "$dir/embeddings.npy" || ! -s "$dir/sample_target_seqid_map.tsv" ]]; then
    echo "[skip] $label: embeddings or map not present yet"
    return 0
  fi
  echo "--- $label ---"
  "$PYTHON_BIN" "$ROOT/.validate_block.py" "$dir" "$VALIDATE_ROWS"
}

echo
echo "=== validation (python=$PYTHON_BIN, rows=$VALIDATE_ROWS) ==="
validate_block "$EMB_DST"  "tb_profiler (all-gene)"
validate_block "$GWAS_DST" "gwas (intergenic)"
rm -f "$ROOT/.validate_block.py"

echo
echo "[OK] root=$ROOT"
if [[ -n "$SHARE_FROM" ]]; then
  echo "[OK] small inputs shared from $SHARE_FROM"
  exit 0
fi
echo "[TODO] upload from the Mac into:"
echo "        $ROOT/metadata/"
echo "        $ROOT/ml_inputs/tbprofiler_binary/"
echo "        $ROOT/ml_inputs/embeddings/gwas/sample_target_seqid_map.tsv   (intergenic-only)"
echo "        $ROOT/$ANTIBIOTIC/results/holdout_results/   (sample_splits_* + sampled_configs_*)"
