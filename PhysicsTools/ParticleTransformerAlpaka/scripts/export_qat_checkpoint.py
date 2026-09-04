#!/usr/bin/env python3
"""Export a QAT ParticleTransformer state_dict to the Alpaka model format.

The quantization-aware checkpoints produced by the modified ParT.py differ from
the ones the original exporter understands in three ways:

  * every quantized layer has an extra ``.wrapped.`` level in its parameter
    names, and carries ``.observer.amax`` / ``.observer.observations`` buffers;
  * attention uses separate ``q_proj`` / ``k_proj`` / ``v_proj`` rather than a
    fused ``in_proj``; the generated graph wants them concatenated;
  * the activation ranges learned during QAT have to travel with the model,
    because the kernel needs them to quantize activations at run time.

Weights are symmetric INT8 per output channel, exactly as before.  What is new
in the emitted header is that every Linear also carries the activation scale of
its input, so an inference kernel can quantize activations instead of
dequantizing weights.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from weaver_checkpoint import load_state_dict  # noqa: E402


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------- #
# checkpoint normalization
# --------------------------------------------------------------------------- #

QKV = ("q_proj", "k_proj", "v_proj")


def normalize(raw: dict) -> tuple[dict, dict]:
    """Return (tensors, input_scales) in the generated graph's naming.

    ``input_scales`` maps a Linear's canonical name to the symmetric per-tensor
    activation scale observed for its input, or None where the layer was left
    in FP32 during training.
    """
    tensors: dict = {}
    amax: dict = {}

    for name, value in raw.items():
        base = name
        if ".observer.amax" in base:
            amax[base[: -len(".observer.amax")]] = float(value)
            continue
        if ".observer.observations" in base:
            if int(value) == 0:
                raise SystemExit(f"{base} was never calibrated; run calibrate_quantization()")
            continue
        base = base.replace(".wrapped.", ".")
        tensors[base] = value

    # Fuse the three attention projections into the in_proj the graph expects.
    scales: dict = {}
    fused: dict = {}
    for key in list(tensors):
        if not key.endswith(".q_proj.weight"):
            continue
        prefix = key[: -len(".q_proj.weight")]
        for suffix in ("weight", "bias"):
            parts = [tensors[f"{prefix}.{p}.{suffix}"] for p in QKV]
            fused[f"{prefix}.in_proj.{suffix}"] = np.concatenate(parts, axis=0)
        # q, k and v observe the same tensor in a particle block, but in a class
        # block the query is the unnormalized class token while the keys and
        # values come from the normalized sequence.  Both scales are kept; the
        # kernel already evaluates the query as a separate projection.
        have = [amax.get(f"{prefix}.{p}") for p in QKV]
        if all(v is not None for v in have):
            # Q, K and V observed their inputs separately: in a particle block
            # all three see the normalised sequence, but in a class block the
            # query sees the raw class token.  Averaging them into one scale is
            # a systematic error, so all three travel with the layer.
            # Stored as (key, query, value): the kernel's fused projection
            # quantizes its input once with the key/value scale, and the class
            # block quantizes the raw class token separately with the query
            # scale.  Key and value always coincide -- they observe the same
            # tensor -- so nothing is lost by sharing one scale between them.
            scales[f"{prefix}.in_proj"] = (have[1] / 127.0, have[0] / 127.0, have[2] / 127.0)
        for p in QKV:
            for suffix in ("weight", "bias"):
                tensors.pop(f"{prefix}.{p}.{suffix}", None)

    tensors.update(fused)
    for key, value in amax.items():
        if not any(key.endswith(f".{p}") for p in QKV):
            scales[key] = value / 127.0
    return tensors, scales


def ordered(tensors: dict) -> list:
    """Deterministic emission order: the graph's own traversal order."""
    order = ["cls_token"]
    order += [f"embed.input_bn.{s}" for s in ("weight", "bias", "running_mean", "running_var")]
    for i in (0, 3, 6):
        order += [f"embed.embed.{i}.weight", f"embed.embed.{i}.bias"]
    for i in (1, 4, 7):
        order += [f"embed.embed.{i}.weight", f"embed.embed.{i}.bias"]
    order += [f"pair_embed.embed.0.{s}" for s in ("weight", "bias", "running_mean", "running_var")]
    for i in (2, 5, 8, 11):
        order += [f"pair_embed.embed.{i}.{s}" for s in ("weight", "bias", "running_mean", "running_var")]
    for i in (1, 4, 7, 10):
        order += [f"pair_embed.embed.{i}.weight", f"pair_embed.embed.{i}.bias"]
    for group, count in (("blocks", 8), ("cls_blocks", 2)):
        for j in range(count):
            n = f"{group}.{j}"
            order += [f"{n}.c_attn", f"{n}.w_resid"]
            for sub in ("pre_attn_norm", "post_attn_norm", "pre_fc_norm", "post_fc_norm"):
                order += [f"{n}.{sub}.weight", f"{n}.{sub}.bias"]
            for sub in ("attn.in_proj", "attn.out_proj", "fc1", "fc2"):
                order += [f"{n}.{sub}.weight", f"{n}.{sub}.bias"]
    order += ["norm.weight", "norm.bias", "fc.0.weight", "fc.0.bias"]
    missing = [k for k in order if k not in tensors]
    if missing:
        raise SystemExit("checkpoint is missing:\n  " + "\n  ".join(missing))
    extra = [k for k in tensors if k not in order and not k.endswith("num_batches_tracked")]
    if extra:
        raise SystemExit("checkpoint has unsupported tensors:\n  " + "\n  ".join(sorted(extra)))
    return order


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("checkpoint", type=pathlib.Path)
    ap.add_argument("--header", type=pathlib.Path, required=True)
    ap.add_argument("--weights", type=pathlib.Path, required=True)
    ap.add_argument("--params", type=pathlib.Path, required=True)
    ap.add_argument("--report", type=pathlib.Path, required=True)
    args = ap.parse_args()

    tensors, scales = normalize(load_state_dict(args.checkpoint))
    order = ordered(tensors)

    weights: list = []
    params: list = []
    locations: dict = {}
    linears: dict = {}
    report = {
        "format": "weaver_part_qat_v1",
        "scheme": "symmetric_int8_per_output_channel_with_activation_scales",
        "checkpoint_sha256": hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
        "tensors": {},
        "activation_scales": {},
    }

    def add_param(name, x):
        x = np.asarray(x, dtype=np.float32).reshape(-1)
        locations[name] = len(params)
        params.extend(x.tolist())

    # Layers that quantization-aware training deliberately left in FP32 must
    # stay in FP32 here too.  The classifier is 256 weights sitting directly
    # under the softmax, and quantizing it shifts the output by far more than
    # every other layer combined -- it was never trained to tolerate it.
    float_weight_layers = {"fc.0"}

    for name in order:
        x = tensors[name]
        if name.endswith(".weight") and name[:-7] in float_weight_layers:
            a = np.asarray(x)
            add_param(name, a.reshape(-1))          # row-major [out][in]
            add_param(name[:-6] + "bias", tensors[name[:-6] + "bias"])
            linears[name[:-7]] = ("float", a.shape[1], a.shape[0],
                                  locations[name], locations[name[:-6] + "bias"])
            report["tensors"][name] = {"shape": list(a.shape), "storage": "float32"}
        elif name.endswith(".weight") and np.asarray(x).ndim >= 2:
            w = np.asarray(x).reshape(np.asarray(x).shape[0], -1).astype(np.float32)
            scale = np.maximum(np.max(np.abs(w), 1), 1e-12) / 127.0
            q = np.clip(np.rint(w / scale[:, None]), -127, 127).astype(np.int8)
            wo = len(weights)
            weights.extend(q.reshape(-1).tolist())
            so = len(params)
            params.extend(scale.tolist())
            bias = tensors.get(name[:-6] + "bias", np.zeros(w.shape[0], np.float32))
            bo = len(params)
            params.extend(np.asarray(bias, np.float32).reshape(-1).tolist())
            layer = name[:-7]
            # The activation scale of this layer's input.  Layers left in FP32
            # during QAT get 0, which the kernel reads as "dequantize weights,
            # do not quantize the activation".
            act = scales.get(layer, 0.0)
            triple = tuple(float(v) for v in act) if isinstance(act, tuple) else (float(act),) * 3
            ao = len(params)
            params.extend(triple)
            linears[layer] = (w.shape[1], w.shape[0], wo, so, bo, ao)
            e = np.abs(w - q.astype(np.float32) * scale[:, None])
            report["tensors"][name] = {
                "shape": list(np.asarray(x).shape),
                "max_abs_error": float(e.max()),
                "mean_abs_error": float(e.mean()),
            }
            report["activation_scales"][layer] = list(triple)
        elif not name.endswith("num_batches_tracked") and name not in\
                {n + ".bias" for n in float_weight_layers}:
            add_param(name, x)

    weights_array = np.ascontiguousarray(weights, dtype=np.int8)
    params_array = np.ascontiguousarray(params, dtype=np.float32)
    weights_bytes = weights_array.tobytes(order="C")
    params_bytes = params_array.tobytes(order="C")

    def lin(n):
        entry = linears[n]
        if entry[0] == "float":
            _, i, o, w, b = entry
            return f"FloatLinear{{{i},{o},m.params+{w},m.params+{b}}}"
        i, o, w, sc, b, a = entry
        return f"Linear{{{i},{o},m.weights+{w},m.params+{sc},m.params+{b},m.params+{a}}}"

    def norm(n):
        size = np.asarray(tensors[n + ".weight"]).size
        return f"Norm{{{size},1e-5f,m.params+{locations[n+'.weight']},m.params+{locations[n+'.bias']}}}"

    def bn(n):
        size = np.asarray(tensors[n + ".weight"]).size
        return (f"BatchNorm{{{size},1e-5f,m.params+{locations[n+'.weight']},m.params+{locations[n+'.bias']},"
                f"m.params+{locations[n+'.running_mean']},m.params+{locations[n+'.running_var']}}}")

    def block(n):
        return (f"BlockData{{{norm(n+'.pre_attn_norm')},{lin(n+'.attn.in_proj')},{lin(n+'.attn.out_proj')},"
                f"{norm(n+'.post_attn_norm')},{norm(n+'.pre_fc_norm')},{lin(n+'.fc1')},"
                f"{norm(n+'.post_fc_norm')},{lin(n+'.fc2')},m.params+{locations[n+'.c_attn']},"
                f"m.params+{locations[n+'.w_resid']}}}")

    out = [
        "// Generated from a QAT Weaver checkpoint; values live in external binary files; do not edit.",
        "#pragma once",
        "#include <cstddef>",
        "#include <cstdint>",
        "namespace generated {",
        f'inline constexpr char checkpoint_sha256[]="{report["checkpoint_sha256"]}";',
        'inline constexpr char storage_precision[]="int8_weights_fp32_params_with_activation_scales";',
        f'inline constexpr char weights_sha256[]="{sha256_bytes(weights_bytes)}";',
        f'inline constexpr char params_sha256[]="{sha256_bytes(params_bytes)}";',
        f"inline constexpr std::size_t weights_size={weights_array.size};",
        f"inline constexpr std::size_t params_size={params_array.size};",
        # ParticleTransformer defaults to use_pre_activation_pair=True, which
        # strips the GELU after the last pair-embedding BatchNorm.
        "inline constexpr bool pair_final_gelu=false;",
        # c_attn is applied as (batch, token, head, head_dim), so the exported
        # head scaling needs no transposition.
        "inline constexpr bool transpose_scaled_heads=false;",
        "struct ModelView { int8_t const* weights; float const* params; };",
        f"ALPAKA_FN_ACC inline BatchNorm inputBN(ModelView m){{return {bn('embed.input_bn')};}}",
        "ALPAKA_FN_ACC inline Norm embedNorm(ModelView m,int i){switch(i){"
        + "".join(f"case {j}:return {norm('embed.embed.'+str(k))};" for j, k in enumerate((0, 3, 6)))
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline Linear embedLinear(ModelView m,int i){switch(i){"
        + "".join(f"case {j}:return {lin('embed.embed.'+str(k))};" for j, k in enumerate((1, 4, 7)))
        + "default:return {};}}",
        f"ALPAKA_FN_ACC inline BatchNorm pairInputBN(ModelView m){{return {bn('pair_embed.embed.0')};}}",
        "ALPAKA_FN_ACC inline BatchNorm pairBN(ModelView m,int i){switch(i){"
        + "".join(f"case {j}:return {bn('pair_embed.embed.'+str(k))};" for j, k in enumerate((2, 5, 8, 11)))
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline Linear pairLinear(ModelView m,int i){switch(i){"
        + "".join(f"case {j}:return {lin('pair_embed.embed.'+str(k))};" for j, k in enumerate((1, 4, 7, 10)))
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline BlockData block(ModelView m,int i){switch(i){"
        + "".join(f"case {j}:return {block('blocks.'+str(j))};" for j in range(8))
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline BlockData clsBlock(ModelView m,int i){switch(i){"
        + "".join(f"case {j}:return {block('cls_blocks.'+str(j))};" for j in range(2))
        + "default:return {};}}",
        f"ALPAKA_FN_ACC inline Norm finalNorm(ModelView m){{return {norm('norm')};}}",
        f"ALPAKA_FN_ACC inline FloatLinear classifier(ModelView m){{return {lin('fc.0')};}}",
        f"ALPAKA_FN_ACC inline float const* clsToken(ModelView m){{return m.params+{locations['cls_token']};}}",
        "} // namespace generated",
        "",
    ]

    for path in (args.header, args.weights, args.params, args.report):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.header.write_text("\n".join(out))
    args.weights.write_bytes(weights_bytes)
    args.params.write_bytes(params_bytes)
    quantized = sum(1 for v in report["activation_scales"].values() if v[0] > 0)
    report.update({
        "weights": {"dtype": "int8", "elements": int(weights_array.size), "sha256": sha256_bytes(weights_bytes)},
        "params": {"dtype": "float32", "elements": int(params_array.size), "sha256": sha256_bytes(params_bytes)},
        "layers_with_activation_scales": quantized,
    })
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "weights_elements": int(weights_array.size),
        "params_elements": int(params_array.size),
        "dense_tensors": len(report["tensors"]),
        "layers_with_activation_scales": quantized,
        "checkpoint_sha256": report["checkpoint_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
