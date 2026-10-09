#!/usr/bin/env bash
# Install the workflow into the layout its scripts expect at run time.
#
# The repository is organised by stage for reading. At run time the scripts call
# one another by fixed locations -- "$MTB/scripts/x.py", "$PIPE/wrappers/y.sh" --
# exactly as they did when the published results were produced, so they are
# copied there rather than rewritten. MANIFEST.tsv says where each file goes.
#
# Usage:
#     source config/paths.env
#     ./install.sh            # copy
#     ./install.sh --dry-run  # show what would be copied
#     ./install.sh --force    # also overwrite installed files that differ
#
# Without --force an existing file that differs is left alone and reported, so
# installing over a tree that is in use cannot silently change a running analysis.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

DRY=0; FORCE=0
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --force)   FORCE=1 ;;
    *) echo "usage: $0 [--dry-run] [--force]" >&2; exit 2 ;;
  esac
done

for v in MTB PSEAER PIPE; do
  [[ -n "${!v:-}" ]] || { echo "[ERROR] $v is unset -- source config/paths.env first" >&2; exit 1; }
done

dest_dir() {
  case "$1" in
    MTB/scripts)    echo "$MTB/scripts" ;;
    PSEAER/scripts) echo "$PSEAER/scripts" ;;
    PIPE/scripts)   echo "$PIPE/scripts" ;;
    PIPE/wrappers)  echo "$PIPE/wrappers" ;;
    PIPE/00_submit) echo "$PIPE/00_submit" ;;
    *) echo "[ERROR] unknown install target $1" >&2; exit 1 ;;
  esac
}

n_new=0; n_same=0; n_over=0; n_skip=0
while IFS=$'\t' read -r repo_path _orig install_to _ _; do
  [[ "$repo_path" == "repo_path" ]] && continue
  IFS=';' read -ra targets <<< "$install_to"
  for t in "${targets[@]}"; do
    d="$(dest_dir "$t")"
    dst="$d/$(basename "$repo_path")"
    if [[ -e "$dst" ]]; then
      if cmp -s "$repo_path" "$dst"; then n_same=$((n_same + 1)); continue; fi
      if [[ "$FORCE" != 1 ]]; then
        echo "[SKIP] differs, not overwritten: $dst"; n_skip=$((n_skip + 1)); continue
      fi
      n_over=$((n_over + 1))
    else
      n_new=$((n_new + 1))
    fi
    if [[ "$DRY" == 1 ]]; then
      echo "[DRY] $repo_path -> $dst"
    else
      mkdir -p "$d"
      cp -p "$repo_path" "$dst"
    fi
  done
done < MANIFEST.tsv

echo "[DONE] new $n_new  unchanged $n_same  overwritten $n_over  skipped $n_skip"
[[ "$n_skip" -eq 0 ]] || echo "        rerun with --force to overwrite the skipped files"
