#!/usr/bin/env python3
"""Export a supported Weaver ParticleTransformer without quantization.

Every checkpoint tensor is preserved as IEEE-754 float32. Matrix weights and
the remaining parameters are written to separate binary blobs; the generated
C++ header contains only sizes, offsets and device-safe accessors.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from weaver_checkpoint import detect_checkpoint_variant, load_state_dict


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path, help="GeneratedModel.h")
    parser.add_argument("--weights", type=pathlib.Path)
    parser.add_argument("--params", type=pathlib.Path)
    parser.add_argument("--report", type=pathlib.Path)
    parser.add_argument(
        "--model-variant",
        choices=("auto", "legacy_mha", "modern_custom"),
        default="auto",
        help="attention implementation; auto detects it from the original QKV key spelling",
    )
    parser.add_argument(
        "--pair-final-activation",
        choices=("auto", "gelu", "none"),
        default="auto",
        help=(
            "activation after the last pair convolution; auto selects gelu for the "
            "verified legacy-MHA model and none for the modern custom-attention model"
        ),
    )
    args = parser.parse_args()

    state = load_state_dict(args.checkpoint)
    detected_variant = detect_checkpoint_variant(args.checkpoint)
    model_variant = detected_variant if args.model_variant == "auto" else args.model_variant
    if args.model_variant != "auto" and args.model_variant != detected_variant:
        raise SystemExit(
            f"--model-variant={args.model_variant} conflicts with checkpoint dialect {detected_variant}"
        )
    pair_final_gelu = (
        model_variant == "legacy_mha"
        if args.pair_final_activation == "auto"
        else args.pair_final_activation == "gelu"
    )
    required = {
        "cls_token",
        "embed.embed.1.weight",
        "pair_embed.embed.10.weight",
        "blocks.0.attn.in_proj.weight",
        "blocks.0.attn.in_proj.bias",
        "blocks.7.fc2.weight",
        "cls_blocks.1.fc2.weight",
        "fc.0.weight",
    }
    if (
        not required <= state.keys()
        or state["embed.embed.1.weight"].shape != (128, 15)
        or state["fc.0.weight"].shape != (2, 128)
    ):
        raise SystemExit("checkpoint is not a supported 15-input, 2-class Weaver ParT architecture")

    weights_path = args.weights or args.output.with_suffix(".weights.fp32.bin")
    params_path = args.params or args.output.with_suffix(".params.fp32.bin")
    report_path = args.report or args.output.with_suffix(".fp32.json")

    weights: list[float] = []
    params: list[float] = []
    locations: dict[str, int] = {}
    linears: dict[str, tuple[int, int, int, int]] = {}
    consumed: set[str] = set()

    report = {
        "format": "weaver_part_v1",
        "precision": "float32_original_no_quantization",
        "model_variant": model_variant,
        "pair_final_activation": "gelu" if pair_final_gelu else "none",
        "embed_final_activation": "gelu",
        "class_query_pre_normalized": False,
        "attention_head_scaling": (
            "after_output_projection_tbdh" if model_variant == "legacy_mha"
            else "after_output_projection_bthd"
        ),
        "checkpoint_sha256": hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
        "tensors": {},
    }

    for name, tensor in state.items():
        if not (name.endswith(".weight") and tensor.ndim >= 2):
            continue
        matrix = np.ascontiguousarray(tensor.reshape(tensor.shape[0], -1), dtype=np.float32)
        weight_offset = len(weights)
        weights.extend(matrix.reshape(-1).tolist())
        prefix = name[:-7]
        bias_name = prefix + ".bias"
        bias = np.ascontiguousarray(
            state.get(bias_name, np.zeros(matrix.shape[0], dtype=np.float32)), dtype=np.float32
        ).reshape(-1)
        bias_offset = len(params)
        params.extend(bias.tolist())
        linears[prefix] = (matrix.shape[1], matrix.shape[0], weight_offset, bias_offset)
        consumed.add(name)
        if bias_name in state:
            consumed.add(bias_name)
        report["tensors"][name] = {
            "shape": list(tensor.shape),
            "dtype": "float32",
            "sha256": sha256_bytes(matrix.tobytes(order="C")),
        }

    for name, tensor in state.items():
        if name in consumed or name.endswith("num_batches_tracked"):
            continue
        value = np.ascontiguousarray(tensor, dtype=np.float32).reshape(-1)
        locations[name] = len(params)
        params.extend(value.tolist())

    weights_array = np.ascontiguousarray(weights, dtype=np.float32)
    params_array = np.ascontiguousarray(params, dtype=np.float32)
    weights_bytes = weights_array.tobytes(order="C")
    params_bytes = params_array.tobytes(order="C")

    def linear(name: str) -> str:
        inputs, outputs, weight_offset, bias_offset = linears[name]
        return f"Linear{{{inputs},{outputs},m.weights+{weight_offset},m.params+{bias_offset}}}"

    def norm(name: str) -> str:
        return (
            f"Norm{{{state[name + '.weight'].size},1e-5f,"
            f"m.params+{locations[name + '.weight']},m.params+{locations[name + '.bias']}}}"
        )

    def batch_norm(name: str) -> str:
        return (
            f"BatchNorm{{{state[name + '.weight'].size},1e-5f,"
            f"m.params+{locations[name + '.weight']},m.params+{locations[name + '.bias']},"
            f"m.params+{locations[name + '.running_mean']},m.params+{locations[name + '.running_var']}}}"
        )

    def block(name: str) -> str:
        return (
            f"BlockData{{{norm(name + '.pre_attn_norm')},{linear(name + '.attn.in_proj')},"
            f"{linear(name + '.attn.out_proj')},{norm(name + '.post_attn_norm')},"
            f"{norm(name + '.pre_fc_norm')},{linear(name + '.fc1')},"
            f"{norm(name + '.post_fc_norm')},{linear(name + '.fc2')},"
            f"m.params+{locations[name + '.c_attn']},m.params+{locations[name + '.w_resid']}}}"
        )

    checkpoint_hash = report["checkpoint_sha256"]
    lines = [
        "// Generated from a Weaver checkpoint without quantization; do not edit.",
        "#pragma once",
        "#include <cstddef>",
        "namespace generated {",
        f'inline constexpr char checkpoint_sha256[]="{checkpoint_hash}";',
        'inline constexpr char storage_precision[]="float32";',
        f'inline constexpr char weights_sha256[]="{sha256_bytes(weights_bytes)}";',
        f'inline constexpr char params_sha256[]="{sha256_bytes(params_bytes)}";',
        f"inline constexpr std::size_t weights_size={weights_array.size};",
        f"inline constexpr std::size_t params_size={params_array.size};",
        f"inline constexpr bool pair_final_gelu={'true' if pair_final_gelu else 'false'};",
        f"inline constexpr bool transpose_scaled_heads={'true' if model_variant == 'legacy_mha' else 'false'};",
        "struct ModelView { float const* weights; float const* params; };",
        f"ALPAKA_FN_ACC inline BatchNorm inputBN(ModelView m){{return {batch_norm('embed.input_bn')};}}",
        "ALPAKA_FN_ACC inline Norm embedNorm(ModelView m,int i){switch(i){"
        + "".join(
            f"case {index}:return {norm('embed.embed.' + str(layer))};"
            for index, layer in enumerate((0, 3, 6))
        )
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline Linear embedLinear(ModelView m,int i){switch(i){"
        + "".join(
            f"case {index}:return {linear('embed.embed.' + str(layer))};"
            for index, layer in enumerate((1, 4, 7))
        )
        + "default:return {};}}",
        f"ALPAKA_FN_ACC inline BatchNorm pairInputBN(ModelView m){{return {batch_norm('pair_embed.embed.0')};}}",
        "ALPAKA_FN_ACC inline BatchNorm pairBN(ModelView m,int i){switch(i){"
        + "".join(
            f"case {index}:return {batch_norm('pair_embed.embed.' + str(layer))};"
            for index, layer in enumerate((2, 5, 8, 11))
        )
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline Linear pairLinear(ModelView m,int i){switch(i){"
        + "".join(
            f"case {index}:return {linear('pair_embed.embed.' + str(layer))};"
            for index, layer in enumerate((1, 4, 7, 10))
        )
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline BlockData block(ModelView m,int i){switch(i){"
        + "".join(f"case {index}:return {block('blocks.' + str(index))};" for index in range(8))
        + "default:return {};}}",
        "ALPAKA_FN_ACC inline BlockData clsBlock(ModelView m,int i){switch(i){"
        + "".join(f"case {index}:return {block('cls_blocks.' + str(index))};" for index in range(2))
        + "default:return {};}}",
        f"ALPAKA_FN_ACC inline Norm finalNorm(ModelView m){{return {norm('norm')};}}",
        f"ALPAKA_FN_ACC inline Linear classifier(ModelView m){{return {linear('fc.0')};}}",
        f"ALPAKA_FN_ACC inline float const* clsToken(ModelView m){{return m.params+{locations['cls_token']};}}",
        "} // namespace generated",
        "",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    weights_path.parent.mkdir(parents=True, exist_ok=True)
    params_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines))
    weights_path.write_bytes(weights_bytes)
    params_path.write_bytes(params_bytes)
    report.update(
        {
            "weights": {
                "path": str(weights_path),
                "elements": int(weights_array.size),
                "bytes": len(weights_bytes),
                "sha256": sha256_bytes(weights_bytes),
            },
            "params": {
                "path": str(params_path),
                "elements": int(params_array.size),
                "bytes": len(params_bytes),
                "sha256": sha256_bytes(params_bytes),
            },
        }
    )
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "header": str(args.output),
                "weights": str(weights_path),
                "params": str(params_path),
                "report": str(report_path),
                "precision": "float32",
                "model_variant": model_variant,
                "pair_final_activation": "gelu" if pair_final_gelu else "none",
                "weights_elements": int(weights_array.size),
                "params_elements": int(params_array.size),
            }
        )
    )


if __name__ == "__main__":
    main()
