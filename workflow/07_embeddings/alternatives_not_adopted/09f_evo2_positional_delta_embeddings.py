#!/usr/bin/env python3
"""09f: evo2 embeddings as per-position deltas from a per-target anchor allele.

Why this exists
---------------
09b mean-pools the residual stream over the whole window. For a 2.2 kb gene whose
resistance signal is a single substitution (katG S315T), that divides the mutated
position's contribution by ~2223. The layer sweep found no better hook layer, so
pooling is the remaining lever.

Subtracting an anchor does NOT by itself help: mean(h_S) - mean(h_A) ==
mean(h_S - h_A), so a mean over deltas is exactly the pooled delta that
09e_build_delta_embeddings.py already produced -- and a per-target constant
offset leaves random-forest splits unchanged, which is why that gave nothing.
Summing rather than averaging is likewise a constant rescale for fixed-length
alleles, so it changes nothing either.

The signal only survives a NONLINEAR reduction over positions. This script emits
three, each as its own variant directory so they can be scored side by side:

  sum      sum_i delta_i          -- the control. Should reproduce 09b's ranking
                                     almost exactly; if it does not, the delta
                                     path is broken and the rest is meaningless.
  maxabs   per dimension, the delta at the position of largest |delta| --
                                     one substitution's deviation, undiluted.
  topk     mean of delta_i over the k positions with largest ||delta_i||_2 --
                                     as above but averaged over a small
                                     neighbourhood, which is more robust when
                                     the model spreads a variant's effect.

Anchor
------
Default is the modal allele of each target (most sample occurrences), not the
H37Rv reference: the modal allele is wild-type in THIS population, so
lineage-fixed differences from H37Rv are not treated as signal. That is the
katG R463L phylogenetic-marker trap recorded in EVALUATION_PROTOCOL.md section 5.

Indels
------
Positions only correspond when lengths match. Alleles whose length differs from
their anchor are reduced over the common prefix (min length) and listed in
length_mismatch.txt with a count printed at the end. If that count is large the
approach needs real alignment and these numbers should not be trusted.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REDUCTIONS = ("sum", "maxabs", "topk")


def load_09b(path: Path):
    """Reuse 09b's model loading, layer discovery and sequence sanitising."""
    spec = importlib.util.spec_from_file_location("evo2_09b", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"[ERROR] cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    sid, buf = None, []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf)
                sid, buf = line[1:].split()[0], []
            elif sid is not None:
                buf.append(line.strip())
    if sid is not None:
        seqs[sid] = "".join(buf)
    return seqs


def target_assignment(map_path: Path) -> tuple[dict[str, str], dict[str, Counter]]:
    """seq_id -> its dominant target, and per-target occurrence counts."""
    counts: dict[str, Counter] = defaultdict(Counter)
    per_seq: dict[str, Counter] = defaultdict(Counter)
    usecols = ["target_id", "seq_id"]
    for chunk in pd.read_csv(map_path, sep="\t", usecols=usecols, dtype=str, chunksize=2_000_000):
        for t, s in zip(chunk["target_id"], chunk["seq_id"]):
            counts[t][s] += 1
            per_seq[s][t] += 1
    multi = sum(1 for s, c in per_seq.items() if len(c) > 1)
    if multi:
        print(f"[WARN] {multi:,} seq_ids appear under more than one target; "
              "each is assigned to its most frequent one")
    assign = {s: c.most_common(1)[0][0] for s, c in per_seq.items()}
    return assign, counts


def hidden_states(mod, model, seq: str, device: str, hook_module) -> np.ndarray:
    """[L, D] float32 residual stream for one sequence."""
    import torch

    token_ids = model.tokenizer.tokenize(seq)
    input_ids = torch.tensor(token_ids, dtype=torch.int64).unsqueeze(0).to(device)
    captured = {}

    def _hook(_m, _i, out):
        captured["t"] = out[0] if isinstance(out, tuple) else out

    handle = hook_module.register_forward_hook(_hook)
    try:
        with torch.no_grad():
            model(input_ids)
    finally:
        handle.remove()
    t = captured.get("t")
    if t is None:
        raise RuntimeError("hook captured nothing")
    if t.dim() == 3:
        t = t[0]
    return t.detach().to(torch.float32).cpu().numpy()


def reduce_delta(delta: np.ndarray, topk: int) -> dict[str, np.ndarray]:
    """delta is [L, D]; return one vector per reduction."""
    out = {"sum": delta.sum(axis=0)}

    idx = np.abs(delta).argmax(axis=0)                      # [D]
    out["maxabs"] = delta[idx, np.arange(delta.shape[1])]   # signed value at that position

    k = min(topk, delta.shape[0])
    mag = np.linalg.norm(delta, axis=1)                     # [L]
    top = np.argpartition(-mag, k - 1)[:k] if k < delta.shape[0] else np.arange(delta.shape[0])
    out["topk"] = delta[top].mean(axis=0)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-fasta", required=True)
    ap.add_argument("--map-tsv", required=True,
                    help="sample_target_seqid_map.tsv, for grouping and anchor choice")
    ap.add_argument("--out-root", required=True,
                    help="one subdirectory per reduction is created here")
    ap.add_argument("--script-09b", default="",
                    help="path to 09b_embed_with_evo2_hidden.py (default: same dir)")
    ap.add_argument("--model-id", default="arcinstitute/evo2_7b")
    ap.add_argument("--device", default="auto")
    ap.add_argument("--hook-layer", default="norm")
    ap.add_argument("--fp8-mode", default="disable")
    ap.add_argument("--topk", type=int, default=8)
    ap.add_argument("--anchor", choices=["modal"], default="modal")
    ap.add_argument("--reductions", default=",".join(REDUCTIONS))
    ap.add_argument("--max-records", type=int, default=0)
    args = ap.parse_args()

    want = [r.strip() for r in args.reductions.split(",") if r.strip()]
    bad = [r for r in want if r not in REDUCTIONS]
    if bad:
        raise SystemExit(f"[ERROR] unknown reduction(s) {bad}; choose from {REDUCTIONS}")

    b09 = Path(args.script_09b) if args.script_09b else \
        Path(__file__).resolve().parent / "09b_embed_with_evo2_hidden.py"
    mod = load_09b(b09)

    seqs = read_fasta(Path(args.input_fasta))
    print(f"[INFO] sequences in fasta: {len(seqs):,}")
    assign, counts = target_assignment(Path(args.map_tsv))

    ids = [s for s in seqs if s in assign]
    missing = len(seqs) - len(ids)
    if missing:
        print(f"[WARN] {missing:,} fasta records absent from the map; skipped")
    ids.sort()
    if args.max_records:
        ids = ids[: args.max_records]
    by_target: dict[str, list[str]] = defaultdict(list)
    for s in ids:
        by_target[assign[s]].append(s)
    print(f"[INFO] targets: {len(by_target):,}   sequences to embed: {len(ids):,}")

    import torch  # noqa: F401
    fp8_mode = mod.normalize_fp8_mode(args.fp8_mode)
    Evo2, label, _ = mod.try_import_evo2()
    device = mod.resolve_device(args.device, torch.cuda.is_available())
    print(f"[INFO] evo2 import: {label}   device={device}")
    model, load_call, _attempts = mod.try_load_model(Evo2, args.model_id, device)
    print(f"[INFO] loaded via: {load_call}")
    # 09b disables transformer_engine's fp8 after loading, not before; the hook
    # output differs if this is skipped, so the delta would not be comparable.
    if fp8_mode == "disable":
        mod.force_disable_te_fp8()
        print("[INFO] transformer_engine fp8 force-disabled")
    layers = mod.discover_layers(model)
    hook_name = args.hook_layer or layers.get("last_block", "")
    hook_module = mod._get_module_by_name(model, hook_name)
    if hook_module is None:
        raise SystemExit(f"[ERROR] could not resolve hook layer {hook_name!r}")
    print(f"[INFO] Hook target: {hook_name} -> {type(hook_module).__name__}")

    dim = None
    vecs: dict[str, dict[str, np.ndarray]] = {r: {} for r in want}
    mismatched: list[str] = []
    done = 0

    for target, members in sorted(by_target.items()):
        occ = counts.get(target, Counter())
        anchor_id = max(members, key=lambda s: (occ.get(s, 0), s))
        anchor_seq, _ = mod.sanitize_sequence(seqs[anchor_id])
        h_anchor = hidden_states(mod, model, anchor_seq, device, hook_module)
        if dim is None:
            dim = int(h_anchor.shape[1])
            print(f"[INFO] hidden dim = {dim}")

        for sid in members:
            seq, _ = mod.sanitize_sequence(seqs[sid])
            if sid == anchor_id:
                h = h_anchor
            else:
                h = hidden_states(mod, model, seq, device, hook_module)
            if h.shape[0] != h_anchor.shape[0]:
                mismatched.append(sid)
                n = min(h.shape[0], h_anchor.shape[0])
                delta = h[:n] - h_anchor[:n]
            else:
                delta = h - h_anchor
            red = reduce_delta(delta, args.topk)
            for r in want:
                vecs[r][sid] = red[r].astype(np.float32)
            done += 1
            if done % 250 == 0:
                print(f"  {done:,}/{len(ids):,}")

    if dim is None:
        raise SystemExit("[ERROR] nothing embedded")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    out_root = Path(args.out_root)
    for r in want:
        run = out_root / f"delta_{r}" / f"run_{stamp}"
        run.mkdir(parents=True, exist_ok=True)
        arr = np.zeros((len(ids), dim), dtype=np.float32)
        for i, sid in enumerate(ids):
            arr[i] = vecs[r][sid]
        np.save(run / "embeddings.npy", arr)
        pd.DataFrame({"row": np.arange(len(ids)), "seq_id": ids}).to_csv(
            run / "embedding_index.tsv", sep="\t", index=False)
        (run / "manifest.json").write_text(json.dumps({
            "reduction": r, "anchor": args.anchor, "topk": args.topk,
            "hook_layer": hook_name, "model_id": args.model_id,
            "n_sequences": len(ids), "dim": dim,
            "n_length_mismatch": len(mismatched),
            "note": "delta from per-target modal allele, reduced over positions",
        }, indent=2))
        link = out_root / f"delta_{r}" / "latest"
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(run)
        print(f"[OK] delta_{r}: {run}/embeddings.npy shape=({len(ids)}, {dim})")

    if mismatched:
        (out_root / "length_mismatch.txt").write_text("\n".join(mismatched) + "\n")
        print(f"[WARN] {len(mismatched):,}/{len(ids):,} alleles differ in length from their "
              f"anchor and were reduced over the common prefix; see length_mismatch.txt")
    else:
        print("[OK] every allele matched its anchor length; no truncation applied")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
