#!/usr/bin/env python3
"""Export a restricted ParticleTransformer as self-contained Alpaka C++.

Weights are symmetric int8, quantized independently for every output row.
Activations, normalization, softmax and accumulation remain float32.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, pathlib, sys
import numpy as np
import torch
from torch import nn


def load_module(path: pathlib.Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def state_dict_from(path):
    obj = torch.load(path, map_location="cpu", weights_only=False)
    sd = obj.get("state_dict", obj) if isinstance(obj, dict) else obj.state_dict()
    return {k.removeprefix("module."): v for k, v in sd.items()}


def validate(m):
    errors = []
    if getattr(m, "for_segmentation", False): errors.append("segmentation output")
    if getattr(m, "pair_embed", None) is not None: errors.append("pair_embed / pairwise attention bias")
    if getattr(m, "include_global_token", False): errors.append("include_global_token")
    for name, x in m.named_modules():
        if isinstance(x, nn.RMSNorm): errors.append(f"RMSNorm at {name}")
        if x.__class__.__name__ == "SwiGLUFFN": errors.append(f"SwiGLU at {name}")
        if x.__class__.__name__ == "Attention":
            if x.q_norm.__class__ is not nn.Identity: errors.append(f"Q/K norm at {name}")
            if x.headwise_attn_output_gate or x.elementwise_attn_output_gate: errors.append(f"attention gate at {name}")
        if x.__class__.__name__ == "Block":
            if x.c_mask is not None or x.c_attn is not None or x.w_resid is not None:
                errors.append(f"learned attention/residual scaling at {name}")
            if x.fc1_g is not None: errors.append(f"SwiGLU block at {name}")
            if x.post_attn_norm.__class__ is not nn.Identity or x.post_fc_norm.__class__ is not nn.Identity:
                errors.append(f"post projection norm at {name}")
    if errors:
        raise SystemExit("Unsupported graph:\n  - " + "\n  - ".join(sorted(set(errors))))


def ident(s):
    return "w_" + "".join(c if c.isalnum() else "_" for c in s)


def floats(a):
    return ",".join(f"{float(x):.9g}f" for x in np.asarray(a).reshape(-1))


def ints(a):
    return ",".join(str(int(x)) for x in np.asarray(a).reshape(-1))


def emit_array(lines, ctype, name, a):
    lines.append(f"inline constexpr {ctype} {name}[{a.size}] = {{{ints(a) if ctype == 'int8_t' else floats(a)}}};")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-file", type=pathlib.Path, required=True)
    ap.add_argument("--factory-file", type=pathlib.Path, required=True)
    ap.add_argument("--factory", default="make_model")
    ap.add_argument("--checkpoint", type=pathlib.Path, required=True)
    ap.add_argument("--output", type=pathlib.Path, required=True)
    ap.add_argument("--max-particles", type=int, required=True)
    ap.add_argument("--report", type=pathlib.Path)
    args = ap.parse_args()
    load_module(args.model_file.resolve(), "ParT")
    factory = load_module(args.factory_file.resolve(), "part_factory")
    model = getattr(factory, args.factory)()
    model.load_state_dict(state_dict_from(args.checkpoint), strict=True)
    model.eval(); validate(model)

    linears = [(n, x) for n, x in model.named_modules() if isinstance(x, nn.Linear)]
    norms = [(n, x) for n, x in model.named_modules() if isinstance(x, nn.LayerNorm)]
    report = {"scheme":"symmetric_int8_per_output_channel_fp32_activation", "tensors":{}}
    out = ["// Generated file; do not edit.", "#pragma once", "#include <cstdint>",
           "namespace generated {"]
    digest = hashlib.sha256(args.checkpoint.read_bytes()).hexdigest()
    out += [f'inline constexpr char checkpoint_sha256[] = "{digest}";',
            f"inline constexpr int max_particles = {args.max_particles};"]
    for n, x in linears:
        w = x.weight.detach().float().numpy()
        scale = np.maximum(np.max(np.abs(w), axis=1), 1e-12) / 127.0
        q = np.clip(np.rint(w / scale[:,None]), -127, 127).astype(np.int8)
        dq = q.astype(np.float32) * scale[:,None]
        base = ident(n)
        emit_array(out, "int8_t", base, q)
        emit_array(out, "float", base+"_scale", scale)
        bias = x.bias.detach().float().numpy() if x.bias is not None else np.zeros(w.shape[0], np.float32)
        emit_array(out, "float", base+"_bias", bias)
        out.append(f"inline constexpr Linear {base}_desc{{{w.shape[1]}, {w.shape[0]}, {base}, {base}_scale, {base}_bias}};")
        e = np.abs(w-dq)
        report["tensors"][n] = {"shape":list(w.shape), "max_abs_error":float(e.max()), "mean_abs_error":float(e.mean())}
    for n, x in norms:
        base=ident(n)
        emit_array(out,"float",base+"_gamma",x.weight.detach().float().numpy())
        emit_array(out,"float",base+"_beta",x.bias.detach().float().numpy())
        out.append(f"inline constexpr LayerNorm {base}_desc{{{x.normalized_shape[0]}, {float(x.eps):.9g}f, {base}_gamma, {base}_beta}};")
    if getattr(model,"cls_token",None) is not None:
        emit_array(out,"float","cls_token",model.cls_token.detach().float().numpy())
    embed_dim=model.norm.normalized_shape[0]
    max_dim=max([embed_dim]+[x.in_features for _,x in linears]+[x.out_features for _,x in linears])
    nclasses=model.fc[-1].out_features if model.fc is not None else embed_dim
    out += [f"inline constexpr int embed_dim = {embed_dim};", f"inline constexpr int scratch_dim = {max_dim};",
            f"inline constexpr int num_classes = {nclasses};"]
    out.append("template <typename A> ALPAKA_FN_ACC inline void infer(A const& acc, float const* input, int n, float* output) {")
    out.append("  float a[(max_particles+1)*scratch_dim] = {}, b[(max_particles+1)*scratch_dim] = {};")
    out.append("  float q[(max_particles+1)*embed_dim] = {}, k[(max_particles+1)*embed_dim] = {}, v[(max_particles+1)*embed_dim] = {};")
    out.append("  if(n<1)n=1; if(n>max_particles)n=max_particles;")
    # embedding graph
    curdim = model.embed.embed[0].normalized_shape[0] if hasattr(model.embed, "embed") and len(model.embed.embed) else embed_dim
    out.append(f"  for(int t=0;t<n;++t) for(int d=0;d<{curdim};++d) a[t*scratch_dim+d]=input[t*{curdim}+d];")
    if not hasattr(model.embed, "embed"):
        raise SystemExit("Unsupported embedding container")
    seq=list(model.embed.embed)
    i=0
    while i<len(seq):
        if not isinstance(seq[i],nn.LayerNorm) or not isinstance(seq[i+1],nn.Linear):
            raise SystemExit("Embedding must be LayerNorm, Linear, activation groups")
        ln_name=next(n for n,x in norms if x is seq[i]); li_name=next(n for n,x in linears if x is seq[i+1])
        out.append(f"  for(int t=0;t<n;++t){{ part::layerNorm(acc,{ident(ln_name)}_desc,&a[t*scratch_dim],&b[t*scratch_dim]); part::linear(acc,{ident(li_name)}_desc,&b[t*scratch_dim],&a[t*scratch_dim]);")
        act=seq[i+2]
        if isinstance(act,nn.GELU): out.append(f"    for(int d=0;d<{seq[i+1].out_features};++d)a[t*scratch_dim+d]=part::gelu(acc,a[t*scratch_dim+d]); }}")
        elif isinstance(act,nn.ReLU): out.append(f"    for(int d=0;d<{seq[i+1].out_features};++d)if(a[t*scratch_dim+d]<0)a[t*scratch_dim+d]=0; }}")
        else: raise SystemExit(f"Unsupported embedding activation {type(act)}")
        curdim=seq[i+1].out_features; i+=3

    def emit_block(block, prefix, cls=False):
        pa=ident(prefix+".pre_attn_norm"); qn=ident(prefix+".attn.q_proj"); kn=ident(prefix+".attn.k_proj"); vn=ident(prefix+".attn.v_proj"); on=ident(prefix+".attn.out_proj")
        pf=ident(prefix+".pre_fc_norm"); f1=ident(prefix+".fc1"); f2=ident(prefix+".fc2")
        if cls:
            out.append("  // class-attention block")
            out.append(f"  for(int d=0;d<embed_dim;++d)b[d]=a[n*scratch_dim+d];")
            out.append(f"  part::layerNorm(acc,{pa}_desc,b,q); part::linear(acc,{qn}_desc,q,q);")
            out.append(f"  part::layerNorm(acc,{pa}_desc,b,k); part::linear(acc,{kn}_desc,k,k); part::linear(acc,{vn}_desc,k,v);")
            out.append(f"  for(int t=0;t<n;++t){{part::layerNorm(acc,{pa}_desc,&a[t*scratch_dim],&b[t*scratch_dim]);part::linear(acc,{kn}_desc,&b[t*scratch_dim],&k[(t+1)*embed_dim]);part::linear(acc,{vn}_desc,&b[t*scratch_dim],&v[(t+1)*embed_dim]);}}")
            out.append(f"  part::attention(acc,q,k,v,1,n+1,{block.num_heads},{on}_desc,b);")
            out.append("  for(int d=0;d<embed_dim;++d)a[n*scratch_dim+d]+=b[d];")
            target="&a[n*scratch_dim]"; count="1"
        else:
            out.append(f"  for(int t=0;t<n;++t)part::layerNorm(acc,{pa}_desc,&a[t*scratch_dim],&b[t*scratch_dim]);")
            out.append(f"  for(int t=0;t<n;++t){{part::linear(acc,{qn}_desc,&b[t*scratch_dim],&q[t*embed_dim]);part::linear(acc,{kn}_desc,&b[t*scratch_dim],&k[t*embed_dim]);part::linear(acc,{vn}_desc,&b[t*scratch_dim],&v[t*embed_dim]);}}")
            out.append(f"  part::attention(acc,q,k,v,n,n,{block.num_heads},{on}_desc,b);")
            out.append("  for(int t=0;t<n;++t)for(int d=0;d<embed_dim;++d)a[t*scratch_dim+d]+=b[t*embed_dim+d];")
            target="a"; count="n"
        out.append(f"  for(int t=0;t<{count};++t){{float* z={target}+t*scratch_dim;part::layerNorm(acc,{pf}_desc,z,b);part::linear(acc,{f1}_desc,b,q);")
        if isinstance(block.act,nn.GELU): out.append(f"    for(int d=0;d<{block.ffn_dim};++d)q[d]=part::gelu(acc,q[d]);")
        elif isinstance(block.act,nn.ReLU): out.append(f"    for(int d=0;d<{block.ffn_dim};++d)if(q[d]<0)q[d]=0;")
        else: raise SystemExit("Unsupported block activation")
        out.append(f"    part::linear(acc,{f2}_desc,q,b);for(int d=0;d<embed_dim;++d)z[d]+=b[d];}}")

    for bi,bk in enumerate(model.blocks): emit_block(bk,f"blocks.{bi}")
    if model.cls_blocks is not None:
        out.append("  for(int d=0;d<embed_dim;++d)a[n*scratch_dim+d]=cls_token[0*embed_dim+d];")
        for bi,bk in enumerate(model.cls_blocks): emit_block(bk,f"cls_blocks.{bi}",True)
        out.append(f"  part::layerNorm(acc,{ident('norm')}_desc,&a[n*scratch_dim],b);")
    else:
        out.append("  for(int d=0;d<embed_dim;++d){b[d]=0;for(int t=0;t<n;++t)b[d]+=a[t*scratch_dim+d];b[d]/=n;} part::layerNorm(acc,w_norm_desc,b,a);")
    if model.fc is None:
        out.append("  for(int d=0;d<embed_dim;++d)output[d]=b[d];")
    else:
        src="b"
        for fi,layer in enumerate(model.fc):
            name=f"fc.{fi}"
            if isinstance(layer,nn.Linear):
                dst="q" if src=="b" else "b"; out.append(f"  part::linear(acc,{ident(name)}_desc,{src},{dst});"); src=dst
            elif isinstance(layer,(nn.GELU,nn.ReLU,nn.Dropout)):
                if isinstance(layer,nn.GELU): out.append(f"  for(int d=0;d<{model.fc[fi-1].out_features};++d){src}[d]=part::gelu(acc,{src}[d]);")
                elif isinstance(layer,nn.ReLU): out.append(f"  for(int d=0;d<{model.fc[fi-1].out_features};++d)if({src}[d]<0){src}[d]=0;")
            else: raise SystemExit(f"Unsupported classifier layer {type(layer)}")
        out.append(f"  for(int d=0;d<num_classes;++d)output[d]={src}[d];")
    out += ["}", "}  // namespace generated", ""]
    # descriptors must precede generated arrays
    prelude = ['#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelTypes.h"']
    out[3:3] = prelude
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(out))
    rp = args.report or args.output.with_suffix(".quantization.json")
    rp.write_text(json.dumps(report, indent=2)+"\n")
    print(f"wrote {args.output} and {rp}; checkpoint sha256={digest}")

if __name__ == "__main__": main()
