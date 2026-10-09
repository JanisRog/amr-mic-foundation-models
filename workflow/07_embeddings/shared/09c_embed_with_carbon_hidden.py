#!/usr/bin/env python
"""Embed DNA windows with a DNA foundation model (Carbon, DNABERT-2, ...).

Mirrors the output contract of 09b_embed_with_evo2_hidden.py so downstream
merge/CV scripts consume Carbon embeddings without modification:
    embeddings.npy        float32 [n_sequences, hidden_dim]
    embedding_index.tsv   row / feature_id / seq_id / length / ...
    sequence_ids.txt
    embed_manifest.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import torch


def write_json(path: Path, obj: dict) -> None:
    with path.open("w") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    sid = None
    chunks: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None and chunks:
                    records.append((sid, "".join(chunks)))
                sid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if sid is not None and chunks:
        records.append((sid, "".join(chunks)))
    return records


def sanitize(seq: str, policy: str) -> Tuple[str, int]:
    """Uppercase, U->T, then apply the ambiguous-base policy.

    'keep' preserves N (matching the evo2 path, so the two feature sets stay
    comparable); 'replace' rewrites N to A for tokenizers that cannot represent
    ambiguous bases. probe_ambiguous() decides which is safe at runtime.
    """
    s = str(seq).upper().replace("U", "T")
    keep = {"A", "C", "G", "T"}
    out = []
    replaced = 0
    for ch in s:
        if ch in keep:
            out.append(ch)
        elif policy == "keep":
            out.append("N")
        else:
            out.append("A")
            replaced += 1
    return "".join(out), replaced


def wrap_dna(seq: str, mode: str) -> str:
    """Carbon uses a hybrid tokenizer: 6-mer for DNA, Qwen3 BPE for English.

    The <dna> tag is what routes input down the 6-mer path, so unwrapped
    sequences are silently tokenized as English text.
    """
    if mode == "none":
        return seq
    if mode == "open":
        return f"<dna>{seq}"
    return f"<dna>{seq}</dna>"


def probe_dna_tag(tok, mode: str) -> dict:
    """Verify <dna> is a real special token rather than being split into pieces."""
    out = {"mode": mode}
    try:
        for tag in ("<dna>", "</dna>"):
            ids = tok(tag, add_special_tokens=False)["input_ids"]
            out[f"{tag}_n_tokens"] = len(ids)
            out[f"{tag}_ids"] = ids[:4]
        probe = wrap_dna("ACGTAC" * 4, mode)
        ids = tok(probe, add_special_tokens=False)["input_ids"]
        out["example_n_tokens"] = len(ids)
        out["example_bp"] = 24
        # 24 bp of DNA should cost ~4 six-mer tokens plus the tag(s).
        out["looks_like_6mer_path"] = len(ids) <= 8
    except Exception as exc:  # noqa: BLE001
        out["error"] = repr(exc)
    return out


def probe_ambiguous(tok) -> dict:
    """Check whether the tokenizer can represent N without falling back to <unk>."""
    probe = "ACGTACGTNNACGTACGT"
    try:
        enc = tok(probe, add_special_tokens=False)["input_ids"]
    except Exception as exc:  # noqa: BLE001
        return {"n_representable": False, "reason": f"tokenize failed: {exc!r}"}
    unk = getattr(tok, "unk_token_id", None)
    n_unk = sum(1 for t in enc if unk is not None and t == unk)
    roundtrip = ""
    try:
        roundtrip = tok.decode(enc)
    except Exception:  # noqa: BLE001
        pass
    return {
        "n_representable": (n_unk == 0) and ("N" in roundtrip.upper() or not roundtrip),
        "n_unk_tokens": int(n_unk),
        "n_tokens": len(enc),
        "roundtrip_has_N": "N" in roundtrip.upper(),
    }


def resolve_device(arg: str) -> str:
    if arg == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    return arg


def resolve_dtype(name: str):
    return {"bf16": torch.bfloat16, "fp16": torch.float16, "fp32": torch.float32}[name]


def coerce(value: str):
    low = value.strip().lower()
    if low in ("true", "false"):
        return low == "true"
    if low in ("none", "null"):
        return None
    for cast in (int, float):
        try:
            return cast(value)
        except ValueError:
            pass
    return value


def load_model(
    model_dir: str,
    revision: str,
    trust_remote_code: bool,
    dtype,
    device: str,
    config_override: Optional[List[str]] = None,
):
    from transformers import AutoConfig, AutoModel, AutoModelForCausalLM, AutoTokenizer

    kwargs = {"trust_remote_code": trust_remote_code}
    if revision:
        kwargs["revision"] = revision

    tok = AutoTokenizer.from_pretrained(model_dir, **kwargs)
    cfg = AutoConfig.from_pretrained(model_dir, **kwargs)

    for item in config_override or []:
        if "=" not in item:
            raise SystemExit(f"[ERROR] --config-override expects key=value, got {item!r}")
        key, raw = item.split("=", 1)
        key = key.strip()
        val = coerce(raw)
        print(f"[INFO] config override: {key}={val!r} (was {getattr(cfg, key, '<unset>')!r})")
        setattr(cfg, key, val)

    mode = "AutoModel"
    try:
        model = AutoModel.from_pretrained(model_dir, config=cfg, torch_dtype=dtype, **kwargs)
    except Exception as exc:  # noqa: BLE001
        print(f"[INFO] AutoModel failed ({exc!r}); falling back to AutoModelForCausalLM")
        model = AutoModelForCausalLM.from_pretrained(
            model_dir, config=cfg, torch_dtype=dtype, **kwargs
        )
        mode = "AutoModelForCausalLM"

    model.eval().to(device)
    return tok, cfg, model, mode


def extract_hidden(out, layer: int):
    """Get token hidden states from any HF-style return.

    Causal LMs (evo2, Carbon) return a ModelOutput; DNABERT-2's custom
    bert_layers returns a plain tuple whose first element is the sequence
    output, so `out.last_hidden_state` does not exist there.
    """
    hs = getattr(out, "hidden_states", None)
    if hs is not None and layer != -1:
        return hs[layer]
    lhs = getattr(out, "last_hidden_state", None)
    if lhs is not None:
        return lhs
    if hs is not None:
        return hs[-1]
    if isinstance(out, (tuple, list)) and len(out) > 0:
        first = out[0]
        if hasattr(first, "dim") and first.dim() == 3:
            return first
    raise RuntimeError(f"could not locate token hidden states in {type(out)}")


def pool(hidden: torch.Tensor, mask: torch.Tensor, how: str) -> torch.Tensor:
    mask = mask.unsqueeze(-1).to(hidden.dtype)
    if how == "mean":
        return (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1e-6)
    if how == "last":
        idx = mask.squeeze(-1).sum(dim=1).long() - 1
        return hidden[torch.arange(hidden.size(0), device=hidden.device), idx.clamp(min=0)]
    raise ValueError(f"unknown pooling: {how}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--revision", default="")
    ap.add_argument(
        "--trust-remote-code",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Carbon's hybrid tokenizer ships as custom code; required per model card",
    )
    ap.add_argument(
        "--add-special-tokens",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="model card wraps DNA in <dna> and passes add_special_tokens=False",
    )
    ap.add_argument(
        "--config-override",
        nargs="*",
        default=[],
        help=(
            "key=value overrides applied to the model config before loading, e.g. "
            "attention_probs_dropout_prob=0.1 to route DNABERT-2 off its Triton "
            "flash-attention kernel (inactive under eval(), so numerics are unchanged)"
        ),
    )
    ap.add_argument(
        "--dna-tag",
        default="both",
        choices=["both", "open", "none"],
        help="<dna>..</dna> routes input to the 6-mer path; 'none' tokenizes as English",
    )
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    ap.add_argument("--dtype", default="bf16", choices=["bf16", "fp16", "fp32"])
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--max-length", type=int, default=0, help="0 = infer from model config")
    ap.add_argument("--pooling", default="mean", choices=["mean", "last"])
    ap.add_argument("--layer", type=int, default=-1, help="hidden_states index; -1 = last")
    ap.add_argument("--ambiguous-policy", default="replace", choices=["replace", "keep"])
    ap.add_argument("--max-records", type=int, default=0, help="0 = all; >0 for smoke tests")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "embed_manifest.json"

    info = {
        "script": "09c_embed_with_carbon_hidden.py",
        "model_name": Path(args.model_dir).name,
        "input_fasta": args.input_fasta,
        "model_dir": args.model_dir,
        "revision": args.revision,
        "pooling": args.pooling,
        "layer": args.layer,
        "dtype": args.dtype,
        "ambiguous_policy": args.ambiguous_policy,
        "model_name": Path(args.model_dir).name,
        "status": "running",
    }
    write_json(manifest_path, info)

    records = read_fasta(Path(args.input_fasta))
    if args.max_records > 0:
        records = records[: args.max_records]
    if not records:
        info["status"] = "fail"
        info["error"] = "no records parsed from input FASTA"
        write_json(manifest_path, info)
        raise SystemExit("[ERROR] no records parsed")

    device = resolve_device(args.device)
    dtype = resolve_dtype(args.dtype)
    if device == "cpu":
        dtype = torch.float32

    tok, cfg, model, load_mode = load_model(
        args.model_dir, args.revision, args.trust_remote_code, dtype, device,
        args.config_override,
    )

    max_len = args.max_length or int(
        getattr(cfg, "max_position_embeddings", 0)
        or getattr(tok, "model_max_length", 0)
        or 4096
    )
    if max_len > 1_000_000:  # sentinel used by some tokenizers
        max_len = 4096

    if tok.pad_token_id is None:
        if tok.eos_token is not None:
            tok.pad_token = tok.eos_token
        else:
            tok.add_special_tokens({"pad_token": "[PAD]"})
            model.resize_token_embeddings(len(tok))

    tagp = probe_dna_tag(tok, args.dna_tag)
    info["dna_tag_probe"] = tagp
    print(f"[INFO] dna-tag probe: {tagp}")
    if args.dna_tag != "none" and not tagp.get("looks_like_6mer_path", False):
        print(
            "[WARN] <dna> tag did not produce a compact 6-mer tokenization; "
            "input may be routed to the English BPE path.",
            flush=True,
        )

    amb = probe_ambiguous(tok)
    info["ambiguous_probe"] = amb
    print(f"[INFO] ambiguous-base probe: {amb}")
    if args.ambiguous_policy == "keep" and not amb.get("n_representable"):
        print(
            "[WARN] --ambiguous-policy keep requested but the tokenizer does not "
            "represent N cleanly; N will map to <unk>. Consider --ambiguous-policy replace.",
            flush=True,
        )

    ids = [sid for sid, _ in records]
    seqs_raw = [s for _, s in records]
    clean = [sanitize(s, args.ambiguous_policy) for s in seqs_raw]
    seqs = [c[0] for c in clean]
    n_replaced = [c[1] for c in clean]
    lengths = [len(s) for s in seqs]

    # Non-overlapping 6-mer tokenizer: a length not divisible by 6 leaves a
    # remainder the tokenizer must handle outside the clean 6-mer path.
    n_unaligned = sum(1 for L in lengths if L % 6 != 0)
    info["n_sequences_not_multiple_of_6"] = int(n_unaligned)
    if n_unaligned:
        print(
            f"[INFO] {n_unaligned}/{len(lengths)} sequences are not a multiple of 6 bp "
            "(trailing partial 6-mer)",
            flush=True,
        )

    # Sort by length so padding waste stays low, then restore input order.
    order = sorted(range(len(seqs)), key=lambda i: lengths[i])
    vecs: dict[int, np.ndarray] = {}
    n_truncated = 0
    tok_lengths: dict[int, int] = {}

    t0 = time.time()
    with torch.inference_mode():
        for start in range(0, len(order), args.batch_size):
            batch_idx = order[start : start + args.batch_size]
            batch = [wrap_dna(seqs[i], args.dna_tag) for i in batch_idx]
            enc = tok(
                batch,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max_len,
                add_special_tokens=args.add_special_tokens,
            )
            enc = {k: v.to(device) for k, v in enc.items()}
            mask = enc.get("attention_mask")
            if mask is None:
                mask = torch.ones_like(enc["input_ids"])

            # hidden_states is (embeddings, layer_1, ..., layer_N), so --layer 16
            # on a 32-layer model is the mid-stack representation.
            want_states = args.layer != -1 or load_mode != "AutoModel"
            if want_states:
                try:
                    out = model(**enc, output_hidden_states=True)
                except TypeError:
                    # custom modeling code that does not accept the kwarg
                    if args.layer != -1:
                        raise RuntimeError(
                            "this model does not support output_hidden_states; "
                            "use --layer -1"
                        )
                    out = model(**enc)
            else:
                out = model(**enc)
            hidden = extract_hidden(out, args.layer)

            pooled = pool(hidden.float(), mask, args.pooling).cpu().numpy()
            for j, i in enumerate(batch_idx):
                vecs[i] = pooled[j]
                tok_lengths[i] = int(mask[j].sum().item())
                if int(mask[j].sum().item()) >= max_len:
                    n_truncated += 1

            done = start + len(batch_idx)
            if done % (args.batch_size * 25) == 0 or done == len(order):
                rate = done / max(time.time() - t0, 1e-6)
                print(f"[INFO] {done}/{len(order)} seqs  {rate:.1f} seq/s", flush=True)

    X = np.stack([vecs[i] for i in range(len(seqs))]).astype(np.float32)

    np.save(outdir / "embeddings.npy", X)
    with (outdir / "embedding_index.tsv").open("w") as handle:
        handle.write("row\tfeature_id\tseq_id\tlength\tn_tokens\tsource_count\tn_replaced\tembed_call\n")
        for i, sid in enumerate(ids):
            handle.write(
                f"{i}\t{sid}\t{sid}\t{lengths[i]}\t{tok_lengths.get(i, 0)}\t0\t{n_replaced[i]}\tcarbon:{args.pooling}:{args.layer}\n"
            )
    with (outdir / "sequence_ids.txt").open("w") as handle:
        for sid in ids:
            handle.write(f"{sid}\n")

    info.update(
        {
            "status": "ok",
            "load_mode": load_mode,
            "resolved_device": device,
            "max_length": max_len,
            "n_sequences": int(X.shape[0]),
            "embedding_dim": int(X.shape[1]),
            "n_sequences_truncated": int(n_truncated),
            "n_sequences_with_replacements": int(sum(1 for x in n_replaced if x > 0)),
            "n_total_base_replacements": int(sum(n_replaced)),
            "elapsed_sec": round(time.time() - t0, 1),
            "output_files": {
                "embeddings_npy": str(outdir / "embeddings.npy"),
                "embedding_index_tsv": str(outdir / "embedding_index.tsv"),
                "sequence_ids_txt": str(outdir / "sequence_ids.txt"),
            },
        }
    )
    write_json(manifest_path, info)

    model_name = Path(args.model_dir).name
    print(f"\n[OK] {model_name} embeddings complete: n={X.shape[0]}, dim={X.shape[1]}")
    if n_truncated:
        print(f"[WARN] {n_truncated} sequences hit max_length={max_len} and were truncated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
