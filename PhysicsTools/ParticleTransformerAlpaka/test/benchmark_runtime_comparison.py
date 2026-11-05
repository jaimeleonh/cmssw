#!/usr/bin/env python3
"""Benchmark CUDA graph inference against the C++ TorchScript engine.

Requires the two bridges built by build_runtime_benchmark.sh and CMSSW's torch.
This measures inference on common device inputs, not cmsRun producers
or the PyTorchAlpaka TensorCollection adapter. JSON contains that scope explicitly.
"""
from __future__ import annotations
import argparse
import ctypes as C
import gc
import hashlib
import json
import pathlib
import statistics
import sys
import time

import torch

PACKAGE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE / 'scripts'))
from nnexport.frontends import torchscript as frontend
from nnexport import allocator, writer
from nnexport.torch_numpy import to_torch
from verify_model import read_printed_jet


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def libraries(directory):
    graph = C.CDLL(str(directory / 'particle_transformer_benchmark_graph.so'))
    ts = C.CDLL(str(directory / 'particle_transformer_benchmark_torch.so'))
    graph.graph_create.argtypes = [C.c_char_p]
    graph.graph_create.restype = C.c_void_p
    graph.graph_destroy.argtypes = [C.c_void_p]
    graph.graph_bytes.argtypes = [C.c_void_p]
    graph.graph_bytes.restype = C.c_size_t
    graph.graph_fits_warp.argtypes = [C.c_void_p]
    graph.graph_error.restype = C.c_char_p
    graph.graph_resources.argtypes = [C.c_int,C.POINTER(C.c_longlong)]
    graph.graph_launch.argtypes = [C.c_void_p]*5 + [C.c_int,C.c_int,C.c_int,C.c_void_p]
    # graph_launch has handle + four data pointers, then jets/layout/SoA/stream.
    ts.ts_create.argtypes = [C.c_char_p,C.c_int]
    ts.ts_create.restype = C.c_void_p
    ts.ts_error.restype = C.c_char_p
    ts.ts_bind.argtypes = [C.c_void_p]*6 + [C.c_int,C.c_int,C.c_int]
    ts.ts_run.argtypes = [C.c_void_p,C.c_int]
    ts.ts_destroy.argtypes = [C.c_void_p]
    return graph, ts


def export(source, output, quantize):
    graph, meta = frontend.build(source, max_tokens=16, input_features=15,
        quantize=quantize, arena_floats=10560, arena_int8=4096)
    placement = allocator.allocate(graph,float_capacity=10560,int8_capacity=4096)
    output.write_bytes(writer.serialise(graph,placement,arena_floats=10560,arena_int8=4096))
    return meta


def timed(call, repeats, inner):
    for _ in range(5): call()
    torch.cuda.synchronize()
    torch.cuda.reset_peak_memory_stats()
    baseline = torch.cuda.memory_allocated()
    free_before,total_device=torch.cuda.mem_get_info()
    samples=[]
    for _ in range(repeats):
        begin,end=torch.cuda.Event(enable_timing=True),torch.cuda.Event(enable_timing=True)
        begin.record()
        for _ in range(inner):call()
        end.record();end.synchronize()
        samples.append(begin.elapsed_time(end)/inner)
    return {'median_ms':statistics.median(samples),'min_ms':min(samples),
        'samples_ms':samples,'torch_baseline_bytes':baseline,
        'torch_peak_allocated_bytes':torch.cuda.max_memory_allocated(),
        'torch_peak_reserved_bytes':torch.cuda.max_memory_reserved(),
        'temporary_peak_bytes':torch.cuda.max_memory_allocated()-baseline,
        'device_used_before_bytes':total_device-free_before,
        'device_used_after_bytes':total_device-torch.cuda.mem_get_info()[0]}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',choices=['deepsets','part'],required=True)
    ap.add_argument('--sizes',default='512,8192,56052')
    ap.add_argument('--batches',default='256,1024,4096,16384,56052')
    ap.add_argument('--repeats',type=int,default=5)
    ap.add_argument('--device',type=int,default=0)
    ap.add_argument('--multiplicity',choices=['full','mixed'],default='full')
    ap.add_argument('--quantize',action='store_true')
    ap.add_argument('--soa',action='store_true',help='Use CMSSW-like component-major tensor strides, with jets contiguous')
    ap.add_argument('--tf32',action='store_true')
    ap.add_argument('--omit-full-batch',action='store_true',help='Only time explicitly requested Torch batch sizes')
    ap.add_argument('--skip-native-timing',action='store_true',help='Run native once for correctness, then time only TorchScript')
    ap.add_argument('--library-dir',type=pathlib.Path,default=pathlib.Path('/tmp'))
    ap.add_argument('--output',type=pathlib.Path,required=True)
    args=ap.parse_args()
    if not torch.cuda.is_available():raise SystemExit('CUDA GPU access is required')
    torch.cuda.set_device(args.device)
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=args.tf32
    torch.backends.cudnn.allow_tf32=args.tf32
    sizes=[int(x) for x in args.sizes.split(',')]
    batches=[int(x) for x in args.batches.split(',')]
    if args.repeats<1 or not sizes or not batches or any(x<1 for x in sizes+batches):
        ap.error('sizes, batches and repeats must all be positive')
    source=PACKAGE.parents[1]/'L1TriggerScouting/JetTagging/data'/(
        'softjet_deepsets_cpu.pt' if args.model=='deepsets' else 'softjet_part_unique_cpu.pt')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    graphpath=args.library_dir / ('benchmark_' + args.model + ('_int8' if args.quantize else '_float') + '.nnm')
    meta=export(source,graphpath,args.quantize)
    graph,ts=libraries(args.library_dir)
    gh=graph.graph_create(str(graphpath).encode())
    if not gh:raise RuntimeError(graph.graph_error().decode())
    torch_model_before=torch.cuda.memory_allocated()
    th=ts.ts_create(str(source).encode(),args.device)
    if not th:raise RuntimeError(ts.ts_error().decode())
    torch_model_bytes=torch.cuda.memory_allocated()-torch_model_before
    native_bytes=graph.graph_bytes(gh)
    result={'measurement_date':'2026-10-15','scope':'CUDA inference engine only: production graph headers vs C++ TorchScript; CMSSW adapters, preprocessing and transfers excluded',
        'model':args.model,'source':str(source),'source_sha256':digest(source),
        'graph':str(graphpath),'graph_sha256':digest(graphpath),'precision':'int8_weight_only' if args.quantize else 'float32',
        'input_layout':'soa' if args.soa else 'contiguous','tf32':args.tf32,'multiplicity':args.multiplicity,'torch':torch.__version__,
        'cuda':torch.version.cuda,'gpu':torch.cuda.get_device_name(args.device),
        'device_total_bytes':torch.cuda.get_device_properties(args.device).total_memory,
        'native_uploaded_bytes':native_bytes,'torch_model_initial_allocated_bytes':torch_model_bytes,'seed':20261015,'repeats':args.repeats,
        'records':[],'conversion':{k:v for k,v in meta.items() if k not in ('layers','tensors','activation_scales')}}
    def checkpoint():args.output.write_text(json.dumps(result,indent=2)+'\n')
    try:
        for n in sizes:
            torch.manual_seed(20261015)
            f0,v0,k0=read_printed_jet(PACKAGE/'test/printed_jet.txt',16,15)
            live=int(k0.sum())
            with torch.inference_mode():
                indices=torch.randint(live,(n,16),device='cuda')
                f=to_torch(torch,f0[:live]).cuda()[indices].contiguous()
                v=to_torch(torch,v0[:live]).cuda()[indices].contiguous()
                masks=torch.ones(n,16,device='cuda')
                if args.multiplicity=='mixed':
                    counts=torch.randint(1,17,(n,1),device='cuda')
                    masks=(torch.arange(16,device='cuda')[None,:]<counts).float()
                    f*=masks[:,:,None];v*=masks[:,:,None]
                points=f[:,:,5:7].contiguous()
                if args.soa:
                    pitch=((n+31)//32)*32
                    def as_soa(value):
                        storage=torch.empty(value.shape[2],16,pitch,device='cuda')
                        view=storage.permute(2,1,0)[:n]
                        view.copy_(value)
                        return view
                    f=as_soa(f);v=as_soa(v);points=as_soa(points)
                    storage=torch.empty(16,pitch,device='cuda')
                    view=storage.T[:n];view.copy_(masks);masks=view
                native_output=torch.empty(n,1,device='cuda')
                torch_output=torch.empty_like(native_output)
                del indices
                if ts.ts_bind(th,f.data_ptr(),v.data_ptr(),masks.data_ptr(),points.data_ptr(),torch_output.data_ptr(),n,args.device,int(args.soa)):
                    raise RuntimeError(ts.ts_error().decode())
                layouts=[0,1] if graph.graph_fits_warp(gh) else [0]
                stream=torch.cuda.current_stream().cuda_stream
                inner=20 if args.model=='deepsets' and n<=8192 else 3 if args.model=='deepsets' else 1
                for layout in layouts:
                    def run_native():
                        if graph.graph_launch(gh,f.data_ptr(),v.data_ptr(),masks.data_ptr(),native_output.data_ptr(),n,layout,int(args.soa),stream):
                            raise RuntimeError(graph.graph_error().decode())
                    torch.cuda.empty_cache()
                    if args.skip_native_timing:
                        run_native();torch.cuda.synchronize();continue
                    measurement=timed(run_native,args.repeats,inner)
                    resources=(C.c_longlong*4)()
                    if graph.graph_resources(layout,resources):raise RuntimeError(graph.graph_error().decode())
                    measurement.update(shared_bytes_per_block=resources[0],local_bytes_per_thread=resources[1],registers_per_thread=resources[2],max_resident_blocks_per_sm=resources[3])
                    measurement.update(engine='graph',jets=n,layout='warp' if layout else 'block',
                        native_model_bytes=native_bytes,jets_per_second=n*1000/measurement['median_ms'])
                    result['records'].append(measurement);checkpoint();print(json.dumps(measurement),flush=True)
                used=set()
                for b in batches+([] if args.omit_full_batch else [n]):
                    b=min(b,n)
                    if b in used:continue
                    used.add(b)
                    def run_torch():
                        if ts.ts_run(th,b):raise RuntimeError(ts.ts_error().decode())
                    torch.cuda.empty_cache()
                    try:
                        measurement=timed(run_torch,args.repeats,inner)
                        torch.cuda.synchronize()
                        error=(native_output-torch_output).abs()
                        finite=bool(torch.isfinite(torch_output).all() and torch.isfinite(native_output).all())
                        measurement.update(engine='torchscript_cpp',jets=n,batch=b,
                            jets_per_second=n*1000/measurement['median_ms'],
                            max_abs_difference=error.max().item(),mean_abs_difference=error.mean().item(),finite=finite)
                        if not finite:raise RuntimeError('Nonfinite output')
                        if args.tf32:
                            candidate=torch_output.clone()
                            torch.backends.cuda.matmul.allow_tf32=False
                            torch.backends.cudnn.allow_tf32=False
                            run_torch();torch.cuda.synchronize()
                            difference=(candidate-torch_output).abs()
                            measurement['tf32_vs_torch_fp32_max_abs']=difference.max().item()
                            measurement['tf32_vs_torch_fp32_mean_abs']=difference.mean().item()
                            torch.backends.cuda.matmul.allow_tf32=True
                            torch.backends.cudnn.allow_tf32=True
                            torch_output.copy_(candidate)
                            del candidate,difference
                    except RuntimeError as exc:
                        measurement={'engine':'torchscript_cpp','jets':n,'batch':b,'error':str(exc)[:1500]}
                        torch.cuda.empty_cache()
                    result['records'].append(measurement);checkpoint();print(json.dumps(measurement),flush=True)
                del f,v,masks,points,native_output,torch_output
                gc.collect();torch.cuda.empty_cache()
    finally:
        torch.cuda.synchronize();ts.ts_destroy(th);graph.graph_destroy(gh);checkpoint()

if __name__=='__main__':main()
