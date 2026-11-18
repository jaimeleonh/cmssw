#!/bin/bash
# Build a self-contained TensorRT EPContext + CUDA-fallback ONNX artifact from an
# ONNX model: quantize it (NVIDIA ModelOpt), build a TensorRT engine with trtexec,
# and embed the engine and the quantized graph into one artifact with
# epcontext_packager.py. This is the artifact format consumed by the TensorRT
# backend of cms::Ort::AlpakaSession (Backend::tensorRT): the embedded engine is
# executed by the TensorRT plugin execution provider, and the embedded ONNX graph
# is the CUDA fallback when TensorRT is unavailable.
#
# Requirements:
#   - $EXPORTER: TensorRT installation prefix providing bin/trtexec and lib/
#   - $PYTHON:   interpreter with the onnx, nvidia-modelopt and tensorrt packages
#
# Usage:
#   EXPORTER=/path/to/TensorRT PYTHON=/path/to/venv/bin/python3 \
#     ./epcontext_exporter.sh model.onnx [output.onnx]
#
# Environment (all optional):
#   QLEVEL                 ModelOpt quantization mode (default: int8)
#   JET_SHAPES             trtexec input shapes; must match the batchSize of the
#                          CMSSW run, since SoftJetIdONNX always feeds full
#                          batches of maxBatchSize jets (default: batch 512, 16
#                          particles, feature-major layout)
#   HARDWARE_ARCHITECTURE  GPU compute capability recorded in the metadata
#                          (default: 8.9)
#   EP_SDK_VERSION         TensorRT version recorded in the metadata (default: 11.4.0)
#
# Quantized model, engine and output default to
#   <model>.<qlevel>.onnx, <model>.<qlevel>.plan, <model>.<qlevel>.fallback.onnx
set -euo pipefail

: "${EXPORTER:?Set EXPORTER to the TensorRT installation prefix (bin/trtexec, lib/)}"
PYTHON="${PYTHON:-python3}"
QLEVEL="${QLEVEL:-int8}"
JET_SHAPES="${JET_SHAPES:-pf_features:512x15x16,pf_vectors:512x4x16,pf_mask:512x1x16}"
HARDWARE_ARCHITECTURE="${HARDWARE_ARCHITECTURE:-8.9}"
EP_SDK_VERSION="${EP_SDK_VERSION:-11.4.0}"
SCRIPTS="$(dirname "$(readlink -f "$0")")"

MODEL="${1:?usage: epcontext_exporter.sh <model.onnx> [output.onnx] (see header for the environment variables)}"
BASE="${MODEL%.onnx}"
QUANTIZED_MODEL="${BASE}.${QLEVEL}.onnx"
ENGINE="${BASE}.${QLEVEL}.plan"
OUTPUT="${2:-${BASE}.${QLEVEL}.fallback.onnx}"
TRTEXEC="${EXPORTER}/bin/trtexec"

[[ -f "$MODEL" ]] || { echo "Missing input model: $MODEL" >&2; exit 1; }
[[ -x "$TRTEXEC" ]] || { echo "Missing executable: $TRTEXEC" >&2; exit 1; }
[[ -f "$SCRIPTS/epcontext_packager.py" ]] || { echo "Missing epcontext_packager.py under $SCRIPTS" >&2; exit 1; }

# TensorRT 11: strongly-typed ONNX; this creates a quantized graph.
# ModelOpt uses random calibration data when calibration_data is not provided.
# Preserve FP32 inputs and outputs for compatibility with CMSSW SoA tensors.
"$PYTHON" - "$MODEL" "$QUANTIZED_MODEL" "$QLEVEL" <<'PY'
import sys

import onnx
from modelopt.onnx.quantization import quantize

input_path, output_path, quantize_mode = sys.argv[1:4]

quantize(
    onnx_path=input_path,
    output_path=output_path,
    quantize_mode=quantize_mode,
    calibration_data=None,
    calibration_method="max",
    calibration_eps=["cpu"],
    high_precision_dtype="fp32",
    direct_io_types=False,
)

converted = onnx.load(output_path)
onnx.save(converted, output_path)

# Verify that ModelOpt preserved the original ONNX input/output interface.
original = onnx.load(input_path)

for kind in ("input", "output"):
    before = [
        (value.name, value.type.tensor_type.elem_type)
        for value in getattr(original.graph, kind)
    ]
    after = [
        (value.name, value.type.tensor_type.elem_type)
        for value in getattr(converted.graph, kind)
    ]

    if before != after:
        raise SystemExit(
            f"{kind} interface changed:\n"
            f"  Original:  {before}\n"
            f"  Converted: {after}"
        )

print(f"Quantized ONNX model: {output_path}")
print("ONNX input/output names and dtypes preserved.")
PY

# TODO: add again --hardwareCompatibilityLevel=ampere+
LD_LIBRARY_PATH="$EXPORTER/lib:${LD_LIBRARY_PATH:-}" \
"$TRTEXEC" \
  --onnx="$QUANTIZED_MODEL" \
  --saveEngine="$ENGINE" \
  --minShapes="$JET_SHAPES" \
  --optShapes="$JET_SHAPES" \
  --maxShapes="$JET_SHAPES" \
  --versionCompatible \
  --excludeLeanRuntime \
  --profilingVerbosity=detailed \
  --skipInference

# The packager embeds the plan and base64-encoded quantized graph in one ONNX.
"$PYTHON" "$SCRIPTS/epcontext_packager.py" \
  --model "$QUANTIZED_MODEL" \
  --engine "$ENGINE" \
  --output "$OUTPUT" \
  --hardware-architecture "$HARDWARE_ARCHITECTURE" \
  --ep-sdk-version "$EP_SDK_VERSION"

printf '\nCreated artifact: %s\n' "$OUTPUT"
printf 'Use Backend::tensorRT (onnxBackend=tensorrt) for primary execution and the CMSSW embedded CUDA fallback.\n'
