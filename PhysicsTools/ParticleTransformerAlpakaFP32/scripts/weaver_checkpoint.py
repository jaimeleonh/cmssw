"""Minimal reader for PyTorch zip state_dict checkpoints, with no torch dependency.

The two supported Weaver ParT implementations serialize the combined QKV
projection differently:

* custom Linear module: ``attn.in_proj.weight`` / ``attn.in_proj.bias``;
* torch MultiheadAttention: ``attn.in_proj_weight`` / ``attn.in_proj_bias``.

Both spellings are returned in the first, canonical form so downstream graph
generation and numerical validation use one unambiguous state-dict contract.
"""
from __future__ import annotations
import collections, io, pickle, zipfile
import numpy as np

class _StorageType:
    def __init__(self, dtype): self.dtype=dtype
class _Storage:
    def __init__(self,key,dtype,size): self.key,self.dtype,self.size=key,dtype,size
class _Tensor:
    def __init__(self,storage,offset,shape,stride):
        self.storage,self.offset,self.shape,self.stride=storage,offset,tuple(shape),tuple(stride)
def _rebuild(storage,offset,shape,stride,*_): return _Tensor(storage,offset,shape,stride)

def canonical_tensor_name(name):
    """Map supported implementation aliases to the generated graph names."""
    name=name.removeprefix("mod.")
    name=name.replace(".attn.in_proj_weight", ".attn.in_proj.weight")
    name=name.replace(".attn.in_proj_bias", ".attn.in_proj.bias")
    return name

def detect_checkpoint_variant(path):
    """Identify the supported Weaver attention implementation.

    The spelling is part of the serialized state_dict and therefore provides a
    reliable discriminator without importing torch or the model source:

    * ``attn.in_proj_weight`` is torch.nn.MultiheadAttention (legacy Weaver);
    * ``attn.in_proj.weight`` is Weaver's newer custom Attention module.

    This distinction matters after the output projection: the legacy Block
    uses ``einsum('tbhd,h->tbdh', ...)`` whereas the newer implementation keeps
    the conventional ``bthd`` channel order.
    """
    with zipfile.ZipFile(path) as z:
        pkl=next(n for n in z.namelist() if n.endswith("data.pkl"))
        meta=_Unpickler(io.BytesIO(z.read(pkl))).load()
    names=tuple(name.removeprefix("mod.") for name in meta)
    if any(name=="blocks.0.attn.in_proj_weight" for name in names):
        return "legacy_mha"
    if any(name=="blocks.0.attn.in_proj.weight" for name in names):
        return "modern_custom"
    raise ValueError("checkpoint does not contain a supported combined-QKV attention block")

class _Unpickler(pickle.Unpickler):
    _dtypes={"FloatStorage":np.float32,"DoubleStorage":np.float64,"HalfStorage":np.float16,
             "LongStorage":np.int64,"IntStorage":np.int32,"BoolStorage":np.bool_}
    def find_class(self,module,name):
        if module=="torch._utils" and name.startswith("_rebuild_tensor"): return _rebuild
        if module=="torch" and name in self._dtypes: return _StorageType(self._dtypes[name])
        if module=="torch" and name=="Size": return tuple
        if module=="collections" and name=="OrderedDict": return collections.OrderedDict
        return super().find_class(module,name)
    def persistent_load(self,pid):
        tag,storage_type,key,_location,size=pid
        if tag!="storage": raise ValueError(f"unknown persistent object {tag}")
        return _Storage(key,storage_type.dtype,size)

def load_state_dict(path):
    with zipfile.ZipFile(path) as z:
        pkl=next(n for n in z.namelist() if n.endswith("data.pkl"))
        root=pkl.rsplit("/",1)[0]
        meta=_Unpickler(io.BytesIO(z.read(pkl))).load()
        result={}
        for name,t in meta.items():
            raw=z.read(f"{root}/data/{t.storage.key}")
            base=np.frombuffer(raw,dtype=t.storage.dtype)
            item=base.dtype.itemsize
            view=np.ndarray(t.shape,dtype=base.dtype,buffer=base,
                            offset=t.offset*item,strides=tuple(s*item for s in t.stride))
            canonical=canonical_tensor_name(name)
            if canonical in result:
                raise ValueError(
                    f"checkpoint contains multiple tensors that map to {canonical!r}")
            result[canonical]=view.copy()
        return result
