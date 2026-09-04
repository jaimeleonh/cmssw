#!/usr/bin/env python3
"""Export a supported Weaver ParT state_dict to Alpaka-compatible C++.

Supports checkpoints produced by both the custom-QKV ParT implementation and
the torch.nn.MultiheadAttention implementation.  Their tensor-name difference
is normalized by :mod:`weaver_checkpoint` before architecture validation.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib,sys
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).parent))
from weaver_checkpoint import detect_checkpoint_variant,load_state_dict

def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument("checkpoint",type=pathlib.Path);p.add_argument("output",type=pathlib.Path);p.add_argument("--weights",type=pathlib.Path,help="output INT8 weight blob");p.add_argument("--params",type=pathlib.Path,help="output FP32 scale/bias/norm blob");p.add_argument("--report",type=pathlib.Path);p.add_argument("--model-variant",choices=("auto","legacy_mha","modern_custom"),default="auto");p.add_argument("--pair-final-activation",choices=("auto","gelu","none"),default="auto");a=p.parse_args()
    s=load_state_dict(a.checkpoint);detected_variant=detect_checkpoint_variant(a.checkpoint);model_variant=detected_variant if a.model_variant=="auto" else a.model_variant
    if a.model_variant!="auto" and a.model_variant!=detected_variant:raise SystemExit(f"--model-variant={a.model_variant} conflicts with checkpoint dialect {detected_variant}")
    pair_final_gelu=model_variant=="legacy_mha" if a.pair_final_activation=="auto" else a.pair_final_activation=="gelu"
    required={"cls_token","pair_embed.embed.10.weight","blocks.7.fc2.weight","cls_blocks.1.fc2.weight","fc.0.weight","blocks.0.attn.in_proj.weight","blocks.0.attn.in_proj.bias"}
    if not required<=s.keys() or s["embed.embed.1.weight"].shape!=(128,15) or s["fc.0.weight"].shape!=(2,128):raise SystemExit("checkpoint is not a supported 15-input, 2-class Weaver ParT architecture")
    weights=[];params=[];locations={};linears={};report={"format":"weaver_part_v1","scheme":"symmetric_int8_per_output_channel_fp32_activation","model_variant":model_variant,"pair_final_activation":"gelu" if pair_final_gelu else "none","embed_final_activation":"gelu","class_query_pre_normalized":False,"attention_head_scaling":"after_output_projection_tbdh" if model_variant=="legacy_mha" else "after_output_projection_bthd","checkpoint_sha256":hashlib.sha256(a.checkpoint.read_bytes()).hexdigest(),"tensors":{}}
    def add_param(name,x):
        x=np.asarray(x,dtype=np.float32).reshape(-1);locations[name]=len(params);params.extend(x.tolist())
    for n,x in s.items():
        if n.endswith(".weight") and x.ndim>=2:
            w=x.reshape(x.shape[0],-1).astype(np.float32);scale=np.maximum(np.max(np.abs(w),1),1e-12)/127;q=np.clip(np.rint(w/scale[:,None]),-127,127).astype(np.int8);wo=len(weights);weights.extend(q.reshape(-1).tolist());so=len(params);params.extend(scale.tolist());bias=s.get(n[:-6]+"bias",np.zeros(w.shape[0],np.float32));bo=len(params);params.extend(np.asarray(bias,np.float32).tolist());linears[n[:-7]]=(w.shape[1],w.shape[0],wo,so,bo);e=np.abs(w-q.astype(np.float32)*scale[:,None]);report["tensors"][n]={"shape":list(x.shape),"max_abs_error":float(e.max()),"mean_abs_error":float(e.mean())}
        elif not n.endswith("num_batches_tracked"):
            add_param(n,x)
    weights_path=a.weights or a.output.with_suffix(".weights.int8.bin");params_path=a.params or a.output.with_suffix(".params.fp32.bin");report_path=a.report or a.output.with_suffix(".quantization.json")
    weights_array=np.ascontiguousarray(weights,dtype=np.int8);params_array=np.ascontiguousarray(params,dtype=np.float32);weights_bytes=weights_array.tobytes(order="C");params_bytes=params_array.tobytes(order="C")
    out=["// Generated from Weaver checkpoint; model values are stored in external binary files; do not edit.","#pragma once","#include <cstddef>","#include <cstdint>","namespace generated {"]
    out += [f'inline constexpr char checkpoint_sha256[]="{report["checkpoint_sha256"]}";',f'inline constexpr char storage_precision[]="int8_weights_fp32_params";',f'inline constexpr char weights_sha256[]="{sha256_bytes(weights_bytes)}";',f'inline constexpr char params_sha256[]="{sha256_bytes(params_bytes)}";',f"inline constexpr std::size_t weights_size={weights_array.size};",f"inline constexpr std::size_t params_size={params_array.size};",f"inline constexpr bool pair_final_gelu={'true' if pair_final_gelu else 'false'};",f"inline constexpr bool transpose_scaled_heads={'true' if model_variant=='legacy_mha' else 'false'};"]
    out += ["struct ModelView { int8_t const* weights; float const* params; };"]
    def lin(n):
        i,o,w,sc,b=linears[n];return f"Linear{{{i},{o},m.weights+{w},m.params+{sc},m.params+{b}}}"
    def norm(n):return f"Norm{{{s[n+'.weight'].size},1e-5f,m.params+{locations[n+'.weight']},m.params+{locations[n+'.bias']}}}"
    def bn(n):return f"BatchNorm{{{s[n+'.weight'].size},1e-5f,m.params+{locations[n+'.weight']},m.params+{locations[n+'.bias']},m.params+{locations[n+'.running_mean']},m.params+{locations[n+'.running_var']}}}"
    def block(n):return f"BlockData{{{norm(n+'.pre_attn_norm')},{lin(n+'.attn.in_proj')},{lin(n+'.attn.out_proj')},{norm(n+'.post_attn_norm')},{norm(n+'.pre_fc_norm')},{lin(n+'.fc1')},{norm(n+'.post_fc_norm')},{lin(n+'.fc2')},m.params+{locations[n+'.c_attn']},m.params+{locations[n+'.w_resid']}}}"
    out += [f"ALPAKA_FN_ACC inline BatchNorm inputBN(ModelView m){{return {bn('embed.input_bn')};}}",
            "ALPAKA_FN_ACC inline Norm embedNorm(ModelView m,int i){switch(i){"+"".join(f"case {j}:return {norm('embed.embed.'+str(k))};" for j,k in enumerate((0,3,6)))+"default:return {};}}",
            "ALPAKA_FN_ACC inline Linear embedLinear(ModelView m,int i){switch(i){"+"".join(f"case {j}:return {lin('embed.embed.'+str(k))};" for j,k in enumerate((1,4,7)))+"default:return {};}}",
            f"ALPAKA_FN_ACC inline BatchNorm pairInputBN(ModelView m){{return {bn('pair_embed.embed.0')};}}",
            "ALPAKA_FN_ACC inline BatchNorm pairBN(ModelView m,int i){switch(i){"+"".join(f"case {j}:return {bn('pair_embed.embed.'+str(k))};" for j,k in enumerate((2,5,8,11)))+"default:return {};}}",
            "ALPAKA_FN_ACC inline Linear pairLinear(ModelView m,int i){switch(i){"+"".join(f"case {j}:return {lin('pair_embed.embed.'+str(k))};" for j,k in enumerate((1,4,7,10)))+"default:return {};}}",
            "ALPAKA_FN_ACC inline BlockData block(ModelView m,int i){switch(i){"+"".join(f"case {j}:return {block('blocks.'+str(j))};" for j in range(8))+"default:return {};}}",
            "ALPAKA_FN_ACC inline BlockData clsBlock(ModelView m,int i){switch(i){"+"".join(f"case {j}:return {block('cls_blocks.'+str(j))};" for j in range(2))+"default:return {};}}",
            f"ALPAKA_FN_ACC inline Norm finalNorm(ModelView m){{return {norm('norm')};}}",f"ALPAKA_FN_ACC inline Linear classifier(ModelView m){{return {lin('fc.0')};}}",f"ALPAKA_FN_ACC inline float const* clsToken(ModelView m){{return m.params+{locations['cls_token']};}}","} // namespace generated",""]
    a.output.parent.mkdir(parents=True,exist_ok=True);weights_path.parent.mkdir(parents=True,exist_ok=True);params_path.parent.mkdir(parents=True,exist_ok=True);report_path.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text("\n".join(out));weights_path.write_bytes(weights_bytes);params_path.write_bytes(params_bytes)
    report.update({"weights":{"path":str(weights_path),"dtype":"int8","elements":int(weights_array.size),"bytes":len(weights_bytes),"sha256":sha256_bytes(weights_bytes)},"params":{"path":str(params_path),"dtype":"float32","elements":int(params_array.size),"bytes":len(params_bytes),"sha256":sha256_bytes(params_bytes)}});report_path.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"header":str(a.output),"weights":str(weights_path),"params":str(params_path),"report":str(report_path),"sha256":report["checkpoint_sha256"],"model_variant":model_variant,"pair_final_activation":"gelu" if pair_final_gelu else "none","quantized_tensors":len(report["tensors"]),"weights_elements":int(weights_array.size),"params_elements":int(params_array.size)}))
if __name__=="__main__":main()
