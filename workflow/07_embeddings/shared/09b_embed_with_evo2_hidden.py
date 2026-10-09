#!/usr/bin/env python3
"""
09b: Extract Evo2 sequence embeddings.

Why NOT mlp.l3 (used in 09):
  The mlp.l3 (down-projection) outputs near-zero values (~1e-14) because the
  MLP gate effectively shuts off for short sequences. The 09 script silently
  saved these degenerate embeddings.

What this script extracts instead:
  Priority 1: Forward hook on last transformer block → residual stream (4096-dim)
               This is the model's full internal representation of the sequence.
  Priority 2: Logits (512-dim) — the LM head output, always has real values.

  The logits vector can also be selected explicitly with --embedding-source logits.

  Includes mandatory quality validation: rejects abs_max < 1e-6.
"""
import argparse
import ctypes
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2))


def read_table_auto(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def try_import_evo2():
    tried = []
    candidates = [
        ("from evo2 import Evo2", "from evo2 import Evo2"),
        ("from evo2.models import Evo2", "from evo2.models import Evo2"),
    ]
    for label, stmt in candidates:
        try:
            ns = {}
            exec(stmt, ns)
            return ns["Evo2"], label, tried
        except Exception as e:
            tried.append({"attempt": label, "error": repr(e)})
    raise ImportError(f"Could not import Evo2 class. attempts={tried}")


def normalize_model_name(model_id: str) -> str:
    return model_id.split("/")[-1] if "/" in model_id else model_id


def infer_expected_hidden_dim(model_id: str) -> int:
    mid = str(model_id).lower()
    if re.search(r"(^|[_-])1b([_-]|$)", mid):
        return 2048
    if re.search(r"(^|[_-])7b([_-]|$)", mid):
        return 4096
    return 0


def resolve_device(device_arg: str, cuda_available: bool) -> str:
    if device_arg == "auto":
        return "cuda:0" if cuda_available else "cpu"
    if device_arg == "cuda":
        return "cuda:0"
    return "cpu"


def try_load_model(Evo2, model_id: str, device: str):
    tried = []
    short_name = normalize_model_name(model_id)
    attempts = [
        (f"Evo2('{short_name}', device='{device}')", lambda: Evo2(short_name, device=device)),
        (f"Evo2(model_name='{short_name}', device='{device}')", lambda: Evo2(model_name=short_name, device=device)),
        (f"Evo2('{short_name}')", lambda: Evo2(short_name)),
        (f"Evo2(model_name='{short_name}')", lambda: Evo2(model_name=short_name)),
    ]
    for label, fn in attempts:
        try:
            return fn(), label, tried
        except Exception as e:
            tried.append({"attempt": label, "error": repr(e)})
    raise RuntimeError(f"Could not load model. attempts={tried}")


def force_disable_te_fp8() -> bool:
    import transformer_engine.pytorch as te  # noqa: F401
    if getattr(te.fp8_autocast, "__name__", "") == "_fp8_force_disabled":
        return False
    original = te.fp8_autocast

    def _fp8_force_disabled(*args, **kwargs):
        kwargs["enabled"] = False
        return original(*args, **kwargs)

    _fp8_force_disabled.__name__ = "_fp8_force_disabled"
    te.fp8_autocast = _fp8_force_disabled
    return True


def normalize_fp8_mode(raw: str) -> str:
    s = str(raw).strip().lower()
    if s in {"disable", "disabled", "off", "false", "0", "no"}:
        return "disable"
    if s in {"auto", "on", "enabled", "true", "1", "yes"}:
        return "auto"
    raise ValueError(f"Invalid --fp8-mode '{raw}'.")


def sanitize_sequence(seq: str) -> Tuple[str, int]:
    s = str(seq).upper().replace("U", "T")
    out = []
    replaced = 0
    for ch in s:
        if ch in {"A", "C", "G", "T", "N"}:
            out.append(ch)
        else:
            out.append("N")
            replaced += 1
    return "".join(out), replaced


# ---------------------------------------------------------------------------
#  Embedding quality checks
# ---------------------------------------------------------------------------

def summarize_embedding_array(X: np.ndarray) -> dict:
    arr = np.asarray(X, dtype=np.float64)
    if arr.ndim != 2 or arr.size == 0:
        return {"rank": int(arr.ndim), "shape": list(arr.shape)}
    row_norms = np.sqrt(np.sum(arr * arr, axis=1))
    return {
        "shape": list(arr.shape),
        "global_mean": float(arr.mean()),
        "global_std": float(arr.std()),
        "global_abs_max": float(np.abs(arr).max()),
        "row_norm_min": float(row_norms.min()),
        "row_norm_median": float(np.median(row_norms)),
        "row_norm_max": float(row_norms.max()),
    }


def vec_is_valid(vec_np: np.ndarray, min_abs_max: float = 1e-6) -> bool:
    """Return True if the vector has real signal, False if degenerate."""
    return float(np.abs(vec_np).max()) > min_abs_max


# ---------------------------------------------------------------------------
#  GPU snapshots
# ---------------------------------------------------------------------------

def _parse_cuda_index(device: str) -> int:
    dev = str(device)
    if dev.startswith("cuda:"):
        try:
            return int(dev.split(":", 1)[1])
        except Exception:
            return 0
    return 0


def _bytes_to_gib(x: int) -> float:
    return float(x) / float(1024**3)


def _gpu_snapshot_torch(device: str, phase: str) -> dict:
    snap = {"phase": phase}
    try:
        import torch  # noqa: F401
        if not torch.cuda.is_available():
            snap["torch_cuda_available"] = False
            return snap
        snap["torch_cuda_available"] = True
        idx = _parse_cuda_index(device)
        snap["cuda_index"] = int(idx)
        props = torch.cuda.get_device_properties(idx)
        snap["gpu_name"] = str(getattr(props, "name", "unknown"))
        total_b = int(getattr(props, "total_memory", 0))
        snap["total_memory_bytes"] = total_b
        snap["total_memory_gib"] = _bytes_to_gib(total_b)
        try:
            torch.cuda.synchronize(idx)
        except Exception:
            pass
        for key, fn in [
            ("memory_allocated_bytes", lambda: int(torch.cuda.memory_allocated(idx))),
            ("memory_reserved_bytes", lambda: int(torch.cuda.memory_reserved(idx))),
            ("max_memory_allocated_bytes", lambda: int(torch.cuda.max_memory_allocated(idx))),
            ("max_memory_reserved_bytes", lambda: int(torch.cuda.max_memory_reserved(idx))),
        ]:
            snap[key] = fn()
            snap[key.replace("_bytes", "_gib")] = _bytes_to_gib(snap[key])
        return snap
    except Exception as e:
        snap["torch_probe_error"] = repr(e)
        return snap


def _gpu_snapshot_nvidia_smi(device: str) -> dict:
    snap = {}
    idx = _parse_cuda_index(device)
    cmd = [
        "nvidia-smi",
        "--query-gpu=index,name,memory.total,memory.used,memory.free,utilization.gpu",
        "--format=csv,noheader,nounits",
    ]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if out.returncode != 0:
            snap["nvidia_smi_error"] = (out.stderr or out.stdout).strip() or f"rc={out.returncode}"
            return snap
        rows = [x.strip() for x in (out.stdout or "").splitlines() if x.strip()]
        if not rows:
            snap["nvidia_smi_error"] = "no rows"
            return snap
        chosen = None
        for row in rows:
            parts = [p.strip() for p in row.split(",")]
            if len(parts) < 6:
                continue
            try:
                row_idx = int(parts[0])
            except Exception:
                continue
            if row_idx == idx:
                chosen = parts
                break
        if chosen is None:
            chosen = [p.strip() for p in rows[0].split(",")]
        if len(chosen) >= 6:
            snap["nvidia_smi_index"] = int(float(chosen[0]))
            snap["nvidia_smi_name"] = str(chosen[1])
            snap["nvidia_smi_memory_total_mib"] = int(float(chosen[2]))
            snap["nvidia_smi_memory_used_mib"] = int(float(chosen[3]))
            snap["nvidia_smi_memory_free_mib"] = int(float(chosen[4]))
            snap["nvidia_smi_utilization_gpu_pct"] = int(float(chosen[5]))
        return snap
    except Exception as e:
        snap["nvidia_smi_error"] = repr(e)
        return snap


def collect_gpu_snapshot(device: str, phase: str) -> dict:
    snap = _gpu_snapshot_torch(device, phase)
    snap.update(_gpu_snapshot_nvidia_smi(device))
    return snap


# ---------------------------------------------------------------------------
#  Model introspection
# ---------------------------------------------------------------------------

def discover_layers(model) -> dict:
    """Discover layer names for hook-based extraction."""
    names = []
    for obj in (model, getattr(model, "model", None)):
        if obj is None:
            continue
        try:
            names.extend([n for n, _ in obj.named_modules() if n])
        except Exception:
            pass

    info = {"all_layer_names_count": len(names)}

    # Find the last transformer block (e.g. blocks.31)
    block_names = [n for n in names if re.match(r"^blocks\.\d+$", n)]
    if block_names:
        block_names.sort(key=lambda x: int(x.split(".")[-1]))
        info["last_block"] = block_names[-1]
        info["n_blocks"] = len(block_names)

    # Find norm layers
    norm_names = [n for n in names if n.endswith("norm")]
    if norm_names:
        info["norm_candidates"] = norm_names

    return info


def _get_module_by_name(model, name: str):
    """Resolve a dotted module name to an actual nn.Module."""
    for obj in (model, getattr(model, "model", None)):
        if obj is None:
            continue
        try:
            mod = obj
            for p in name.split("."):
                mod = getattr(mod, p)
            return mod
        except (AttributeError, TypeError):
            continue
    return None


# ---------------------------------------------------------------------------
#  Core extraction
# ---------------------------------------------------------------------------

def _pool_to_np(tensor) -> Optional[np.ndarray]:
    """Mean-pool [batch, seq_len, dim] → [dim] numpy array."""
    import torch  # noqa: F401
    if tensor is None or not torch.is_tensor(tensor):
        return None
    v = tensor.detach()
    if v.ndim >= 3:
        v = v.mean(dim=1)
    elif v.ndim == 1:
        v = v.unsqueeze(0)
    return v.squeeze(0).float().cpu().numpy().astype(np.float32)


def forward_embed_once(
    model, seq: str, device: str,
    hook_module=None, hook_layer_name: str = "",
    expected_dim: int = 0,
    embedding_source: str = "residual",
):
    """
    Extract a sequence embedding.

    embedding_source="residual" uses a forward hook on the configured block and
    falls back to logits if the hook fails. embedding_source="logits" skips hooks
    and intentionally mean-pools the final vocab-size logits vector.
    """
    import torch  # noqa: F401

    token_ids = model.tokenizer.tokenize(seq)
    input_ids = torch.tensor(token_ids, dtype=torch.int64).unsqueeze(0).to(device)

    if embedding_source == "logits":
        with torch.no_grad():
            out_full = model(input_ids)

        logits = out_full
        if isinstance(logits, tuple):
            logits = logits[0]
        if isinstance(logits, (list, tuple)):
            logits = logits[0]
        vec_np = _pool_to_np(logits)
        if vec_np is not None and vec_is_valid(vec_np):
            dim_note = (
                f" dim={vec_np.shape[-1]}"
                if expected_dim <= 0 or int(vec_np.shape[-1]) == int(expected_dim)
                else f" dim={vec_np.shape[-1]} (expected {expected_dim})"
            )
            return vec_np, f"logits mean-pool{dim_note}"
        raise RuntimeError("Logits extraction failed or produced degenerate output.")

    # --- Strategy 1: Forward hook on last block (residual stream, 4096-dim) ---
    if hook_module is not None:
        captured = {}

        def _hook(_mod, _inp, out):
            if isinstance(out, tuple):
                captured["tensor"] = out[0]
            elif torch.is_tensor(out):
                captured["tensor"] = out

        handle = hook_module.register_forward_hook(_hook)
        try:
            with torch.no_grad():
                out_full = model(input_ids)
        finally:
            handle.remove()

        if "tensor" in captured and captured["tensor"] is not None:
            vec_np = _pool_to_np(captured["tensor"])
            if vec_np is not None and vec_is_valid(vec_np):
                dim_ok = (expected_dim <= 0 or int(vec_np.shape[-1]) == int(expected_dim))
                if dim_ok:
                    return vec_np, f"hook({hook_layer_name}) residual-stream dim={vec_np.shape[-1]}"
                # Wrong dim but valid — still report
                return vec_np, f"hook({hook_layer_name}) dim={vec_np.shape[-1]} (expected {expected_dim})"
    else:
        # Need a forward pass for the logits fallback
        with torch.no_grad():
            out_full = model(input_ids)

    # --- Strategy 2: Logits fallback (512-dim, always has real values) ---
    logits = out_full
    if isinstance(logits, tuple):
        logits = logits[0]
    if isinstance(logits, (list, tuple)):
        logits = logits[0]
    vec_np = _pool_to_np(logits)
    if vec_np is not None and vec_is_valid(vec_np):
        return vec_np, f"logits mean-pool dim={vec_np.shape[-1]}"

    raise RuntimeError(
        "All extraction strategies failed or produced degenerate output. "
        f"hook_layer={hook_layer_name}, expected_dim={expected_dim}"
    )


def embed_with_fp8_fallback(model, seq, device, hook_module, hook_layer_name, expected_dim, state):
    try:
        return forward_embed_once(
            model, seq, device, hook_module, hook_layer_name, expected_dim,
            state.get("embedding_source", "residual"),
        )
    except AssertionError as e:
        msg = str(e)
        if ("FP8" in msg or "compute capability" in msg) and not state.get("fp8_disabled", False):
            force_disable_te_fp8()
            state["fp8_disabled"] = True
            state["fp8_disabled_by"] = "auto_fallback_assertion"
            return forward_embed_once(
                model, seq, device, hook_module, hook_layer_name, expected_dim,
                state.get("embedding_source", "residual"),
            )
        raise


def load_records(
    input_fasta: Optional[Path],
    input_tsv: Optional[Path],
    seq_id_col: str,
    seq_col: str,
    count_col: str,
) -> Tuple[List[Tuple[str, str, int]], str]:
    records = []

    if input_tsv is not None:
        if (not input_tsv.exists()) or (not input_tsv.is_file()):
            raise RuntimeError(f"Input TSV does not exist or is not a file: {input_tsv}")
        df = read_table_auto(input_tsv)
        if df.empty:
            return [], "tsv"
        if seq_id_col not in df.columns:
            raise RuntimeError(f"Missing seq id column '{seq_id_col}' in {input_tsv}")
        if seq_col not in df.columns:
            raise RuntimeError(f"Missing sequence column '{seq_col}' in {input_tsv}")

        for _, r in df.iterrows():
            sid = str(r[seq_id_col]).strip()
            seq = str(r[seq_col]).strip()
            if not sid or sid.lower() == "nan":
                continue
            if not seq or seq.lower() == "nan":
                continue
            raw_count = r[count_col] if count_col in df.columns else ""
            try:
                count = int(float(str(raw_count))) if str(raw_count).strip() else 0
            except Exception:
                count = 0
            records.append((sid, seq, count))
        return records, "tsv"

    if input_fasta is None:
        return [], "none"
    if (not input_fasta.exists()) or (not input_fasta.is_file()):
        raise RuntimeError(f"Input FASTA does not exist or is not a file: {input_fasta}")

    try:
        from Bio import SeqIO  # noqa: F401
    except Exception as e:
        raise RuntimeError(f"Biopython is required for FASTA input: {repr(e)}")

    for rec in SeqIO.parse(str(input_fasta), "fasta"):
        sid = str(rec.id).strip()
        seq = str(rec.seq).strip()
        if not sid or not seq:
            continue
        records.append((sid, seq, 0))
    return records, "fasta"


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract Evo2 embeddings from the residual stream (last block output, 4096-dim) "
            "with logits (512-dim) as fallback. Unlike 09 which uses mlp.l3 (near-zero output), "
            "this targets the full residual stream which contains real signal."
        ),
    )
    ap.add_argument("--input-fasta", default="")
    ap.add_argument("--input-tsv", default="")
    ap.add_argument("--seq-id-col", default="seq_id")
    ap.add_argument("--seq-col", default="sequence")
    ap.add_argument("--count-col", default="count")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--model-id", default="arcinstitute/evo2_7b")
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    ap.add_argument("--max-records", type=int, default=0, help="0 = all")
    ap.add_argument(
        "--hook-layer", default="",
        help="Explicit layer name for forward hook (e.g. 'blocks.31'). "
             "Default: auto-detect last block.",
    )
    ap.add_argument(
        "--embedding-source",
        choices=["residual", "logits"],
        default="residual",
        help=(
            "residual = hook the last block and fall back to logits; "
            "logits = intentionally mean-pool the final vocab-size logits vector."
        ),
    )
    ap.add_argument(
        "--expected-embedding-dim", type=int, default=0,
        help="Expected dim. 0 = infer hidden dim for residual or 512 for logits.",
    )
    ap.add_argument(
        "--fp8-mode",
        default=os.environ.get("EVO2_FP8_MODE", "disable"),
    )
    ap.add_argument("--strict-seq", action="store_true")
    args = ap.parse_args()
    fp8_mode = normalize_fp8_mode(args.fp8_mode)

    if bool(args.input_fasta) == bool(args.input_tsv):
        raise RuntimeError("Provide exactly one of --input-fasta or --input-tsv.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "09b_embed_with_evo2_hidden.manifest.json"

    info = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": "09b_embed_with_evo2_hidden.py",
        "description": "Residual-stream embeddings (hook on last block), logits fallback",
        "input_fasta": str(Path(args.input_fasta).resolve()) if args.input_fasta else "",
        "input_tsv": str(Path(args.input_tsv).resolve()) if args.input_tsv else "",
        "seq_id_col": args.seq_id_col,
        "seq_col": args.seq_col,
        "count_col": args.count_col,
        "model_id": args.model_id,
        "device_arg": args.device,
        "max_records": args.max_records,
        "hook_layer_arg": args.hook_layer,
        "embedding_source_arg": args.embedding_source,
        "expected_embedding_dim_arg": int(args.expected_embedding_dim),
        "fp8_mode_resolved": fp8_mode,
        "strict_seq": bool(args.strict_seq),
        "status": "started",
    }

    # --- Environment check ---
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
        import tokenizers  # noqa: F401
        import huggingface_hub  # noqa: F401
        import transformer_engine  # noqa: F401
        import transformer_engine.pytorch  # noqa: F401
        import flash_attn_2_cuda  # noqa: F401
        import vortex  # noqa: F401

        info["python"] = os.path.realpath(os.sys.executable)
        info["torch_version"] = torch.__version__
        info["torch_cuda_available"] = bool(torch.cuda.is_available())
        info["torch_cuda_device_count"] = int(torch.cuda.device_count())
        info["cuda_visible_devices"] = os.environ.get("CUDA_VISIBLE_DEVICES")
        info["versions"] = {
            "transformers": transformers.__version__,
            "tokenizers": tokenizers.__version__,
            "huggingface_hub": huggingface_hub.__version__,
            "transformer_engine": getattr(transformer_engine, "__version__", "unknown"),
        }

        if os.environ.get("CUDART_SO"):
            try:
                ctypes.CDLL(os.environ["CUDART_SO"], mode=ctypes.RTLD_GLOBAL)
                info["cudart_preload"] = "ok"
            except Exception as e:
                info["cudart_preload"] = f"fail: {repr(e)}"

        resolved_device = resolve_device(args.device, info["torch_cuda_available"])
        info["resolved_device"] = resolved_device
        if not resolved_device.startswith("cuda"):
            raise RuntimeError("CUDA is required for Evo2 embedding.")
        info["gpu_snapshot_at_start"] = collect_gpu_snapshot(resolved_device, "env_ready")

    except Exception as e:
        info["status"] = "fail"
        info["error"] = f"env check failed: {repr(e)}"
        write_json(manifest_path, info)
        raise

    # --- Import & load model ---
    try:
        Evo2, import_label, import_attempts = try_import_evo2()
        info["evo2_import"] = {"selected": import_label, "failed_attempts": import_attempts}
    except Exception as e:
        info["status"] = "fail"
        info["error"] = f"Evo2 import failed: {repr(e)}"
        write_json(manifest_path, info)
        raise

    try:
        model, load_call, load_attempts = try_load_model(Evo2, args.model_id, info["resolved_device"])
        info["model_load_call"] = load_call
        info["model_load_failed_attempts"] = load_attempts
        info["gpu_snapshot_after_model_load"] = collect_gpu_snapshot(info["resolved_device"], "after_model_load")
        try:
            import torch  # noqa: F811
            torch.cuda.reset_peak_memory_stats(_parse_cuda_index(info["resolved_device"]))
        except Exception as e:
            info["gpu_peak_reset_error"] = repr(e)
    except Exception as e:
        info["status"] = "fail"
        info["error"] = f"Model load failed: {repr(e)}"
        write_json(manifest_path, info)
        raise

    # --- Discover layers & set up hook ---
    layer_info = discover_layers(model)
    info["discovered_layers"] = layer_info
    print(f"[INFO] Discovered layers: {json.dumps(layer_info, indent=2)}")

    embedding_source = str(args.embedding_source).strip().lower()
    info["embedding_source_resolved"] = embedding_source

    hook_layer_name = args.hook_layer.strip()
    if embedding_source == "logits":
        hook_layer_name = ""
    elif not hook_layer_name:
        hook_layer_name = layer_info.get("last_block", "")
    info["hook_layer_used"] = hook_layer_name

    hook_module = None
    if embedding_source == "logits":
        print("[INFO] Extraction source: logits mean-pool (512-dim vocab vector)")
    elif hook_layer_name:
        hook_module = _get_module_by_name(model, hook_layer_name)
        if hook_module is None:
            print(f"[WARN] Could not resolve hook layer '{hook_layer_name}', will use logits fallback")
        else:
            print(f"[INFO] Hook target: {hook_layer_name} -> {type(hook_module).__name__}")

    expected_dim = int(args.expected_embedding_dim)
    if expected_dim <= 0:
        expected_dim = 512 if embedding_source == "logits" else infer_expected_hidden_dim(args.model_id)
    info["expected_embedding_dim_resolved"] = int(expected_dim)

    # --- Load records ---
    input_fasta = Path(args.input_fasta) if args.input_fasta else None
    input_tsv = Path(args.input_tsv) if args.input_tsv else None
    records, records_from = load_records(
        input_fasta=input_fasta,
        input_tsv=input_tsv,
        seq_id_col=args.seq_id_col,
        seq_col=args.seq_col,
        count_col=args.count_col,
    )
    if args.max_records > 0:
        records = records[: args.max_records]
    if not records:
        info["status"] = "fail"
        info["error"] = "No sequences found in input."
        write_json(manifest_path, info)
        raise RuntimeError(info["error"])

    info["records_from"] = records_from
    info["n_records"] = len(records)

    # --- FP8 ---
    fp8_state = {"fp8_disabled": False, "fp8_disabled_by": "", "embedding_source": embedding_source}
    if fp8_mode == "disable":
        force_disable_te_fp8()
        fp8_state["fp8_disabled"] = True
        fp8_state["fp8_disabled_by"] = "configured_disable_mode"

    # --- Embedding loop ---
    ids = []
    lengths = []
    source_counts = []
    n_replaced_list = []
    calls = []
    vectors = []
    seen_ids = set()

    try:
        for idx, (seq_id, seq_raw, source_count) in enumerate(records):
            if seq_id in seen_ids:
                raise RuntimeError(f"Duplicate sequence id found: {seq_id}")
            seen_ids.add(seq_id)

            seq, n_replaced = sanitize_sequence(seq_raw)
            if args.strict_seq and n_replaced > 0:
                bad_chars = sorted(set(re.sub(r"[ACGTNUacgtnu]", "", str(seq_raw))))
                raise RuntimeError(f"Invalid sequence characters for id={seq_id}: {bad_chars}")
            if not seq:
                raise RuntimeError(f"Empty sequence for id={seq_id}")

            if fp8_mode == "auto":
                vec_np, call = embed_with_fp8_fallback(
                    model, seq, info["resolved_device"],
                    hook_module, hook_layer_name, expected_dim, fp8_state,
                )
            else:
                vec_np, call = forward_embed_once(
                    model, seq, info["resolved_device"],
                    hook_module, hook_layer_name, expected_dim,
                    embedding_source,
                )
            call = f"{call} [fp8_mode={fp8_mode}]"

            if "gpu_snapshot_after_first_forward" not in info:
                info["gpu_snapshot_after_first_forward"] = collect_gpu_snapshot(
                    info["resolved_device"], "after_first_forward"
                )

            if idx == 0:
                abs_max_first = float(np.abs(vec_np).max())
                info["first_embedding_abs_max"] = abs_max_first
                info["first_embedding_call"] = call
                info["first_embedding_dim"] = int(vec_np.shape[0])
                print(f"[INFO] First embedding: dim={vec_np.shape[0]}, "
                      f"abs_max={abs_max_first:.4f}, call={call}")
                if abs_max_first < 1e-6:
                    raise RuntimeError(
                        f"First embedding is degenerate (abs_max={abs_max_first:.3e}). "
                        "Extraction failed to produce real values."
                    )

            ids.append(seq_id)
            lengths.append(len(seq))
            source_counts.append(int(source_count))
            n_replaced_list.append(int(n_replaced))
            calls.append(call)
            vectors.append(vec_np)

            if (idx + 1) % 500 == 0:
                print(f"[INFO] Embedded {idx + 1}/{len(records)} sequences")

        X = np.stack(vectors, axis=0)
        if X.ndim != 2:
            raise RuntimeError(f"Unexpected embedding array rank: shape={X.shape}")

        # Quality validation on full array
        quality = summarize_embedding_array(X)
        info["embedding_quality_summary"] = quality
        abs_max = quality.get("global_abs_max", 0.0)
        if abs_max < 1e-6:
            raise RuntimeError(
                f"DEGENERATE EMBEDDINGS: global_abs_max={abs_max:.3e}. "
                "The extraction produced near-zero values — same bug as mlp.l3. Aborting."
            )

        print(f"[OK] Embedding quality: abs_max={abs_max:.4f}, "
              f"std={quality['global_std']:.4f}, "
              f"row_norm_median={quality['row_norm_median']:.4f}")

    except Exception as e:
        info["status"] = "fail"
        info["fp8_forced_disabled"] = bool(fp8_state["fp8_disabled"])
        info["fp8_disabled_by"] = fp8_state.get("fp8_disabled_by", "")
        info["error"] = f"Embedding loop failed: {repr(e)}"
        write_json(manifest_path, info)
        raise

    # --- Save outputs ---
    np.save(outdir / "embeddings.npy", X)

    with (outdir / "embedding_index.tsv").open("w") as f:
        f.write("row\tfeature_id\tseq_id\tlength\tsource_count\tn_replaced\tembed_call\n")
        for i, (sid, L, cnt, nrep, c) in enumerate(
            zip(ids, lengths, source_counts, n_replaced_list, calls)
        ):
            f.write(f"{i}\t{sid}\t{sid}\t{L}\t{cnt}\t{nrep}\t{c}\n")

    with (outdir / "sequence_ids.txt").open("w") as f:
        for sid in ids:
            f.write(f"{sid}\n")

    # --- Manifest ---
    info["status"] = "ok"
    info["fp8_forced_disabled"] = bool(fp8_state["fp8_disabled"])
    info["fp8_disabled_by"] = fp8_state.get("fp8_disabled_by", "")
    info["n_sequences"] = int(X.shape[0])
    info["embedding_dim"] = int(X.shape[1])

    strategy_counts = {}
    for c in calls:
        key = c.split(" [fp8_mode")[0]
        strategy_counts[key] = strategy_counts.get(key, 0) + 1
    info["extraction_strategy_counts"] = strategy_counts

    info["n_sequences_with_replacements"] = int(sum(x > 0 for x in n_replaced_list))
    info["n_total_base_replacements"] = int(sum(n_replaced_list))
    info["output_files"] = {
        "embeddings_npy": str(outdir / "embeddings.npy"),
        "embedding_index_tsv": str(outdir / "embedding_index.tsv"),
        "sequence_ids_txt": str(outdir / "sequence_ids.txt"),
    }
    info["gpu_snapshot_after_embeddings"] = collect_gpu_snapshot(info["resolved_device"], "after_embeddings")
    write_json(manifest_path, info)

    print(f"\n[OK] Embeddings complete: n={X.shape[0]}, dim={X.shape[1]}")
    print(f"[OK] abs_max={abs_max:.4f}, std={quality['global_std']:.4f}")
    print(f"[OK] Strategies used: {strategy_counts}")
    print(f"[OK] Saved to: {outdir}")
    print(f"[OK] Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
