"""NumPy reference for the exact b_kinadd Weaver ParticleTransformer graph."""
from __future__ import annotations
import math
import numpy as np

def linear(x,w,b): return x@w.T+b
def ln(x,w,b,eps=1e-5):
    return (x-x.mean(-1,keepdims=True))/np.sqrt(x.var(-1,keepdims=True)+eps)*w+b
def bn(x,w,b,mean,var,eps=1e-5):
    shape=(1,)*(x.ndim-1)+(-1,)
    return (x-mean.reshape(shape))/np.sqrt(var.reshape(shape)+eps)*w.reshape(shape)+b.reshape(shape)
def gelu(x):
    erf=np.vectorize(math.erf,otypes=[np.float32])
    return .5*x*(1+erf(x/math.sqrt(2.0)))
def softmax(x):
    e=np.exp(x-np.max(x,axis=-1,keepdims=True)); return e/e.sum(axis=-1,keepdims=True)

def qlinear(x,w,b):
    s=np.maximum(np.max(np.abs(w),axis=1),1e-12)/127
    q=np.clip(np.rint(w/s[:,None]),-127,127).astype(np.int8)
    return x@(q.astype(np.float32)*s[:,None]).T+b

class WeaverParT:
    def __init__(self,state,quantized=False,model_variant="legacy_mha",pair_final_gelu=None):
        if model_variant not in ("legacy_mha","modern_custom"):
            raise ValueError(f"unsupported model variant: {model_variant}")
        self.s=state
        self.lin=qlinear if quantized else linear
        self.model_variant=model_variant
        self.pair_final_gelu=(model_variant=="legacy_mha") if pair_final_gelu is None else pair_final_gelu
    def _norm(self,x,n): return ln(x,self.s[n+".weight"],self.s[n+".bias"])
    def _linear(self,x,n): return self.lin(x,self.s[n+".weight"].reshape(self.s[n+".weight"].shape[0],-1),self.s[n+".bias"])
    def _batchnorm(self,x,n): return bn(x,self.s[n+".weight"],self.s[n+".bias"],self.s[n+".running_mean"],self.s[n+".running_var"])
    def pair(self,v):
        px,py,pz,e=[v[:,i] for i in range(4)]
        pt=np.sqrt(px*px+py*py); rap=.5*np.log(1+2*pz/np.maximum(e-pz,1e-20)); phi=np.arctan2(py,px)
        drap=rap[:,:,None]-rap[:,None,:]; dphi=(phi[:,:,None]-phi[:,None,:]+np.pi)%(2*np.pi)-np.pi
        delta=np.sqrt(drap*drap+dphi*dphi); ptmin=np.minimum(pt[:,:,None],pt[:,None,:])
        f=[np.log(np.maximum(ptmin*delta,1e-8)),np.log(np.maximum(ptmin/np.maximum(pt[:,:,None]+pt[:,None,:],1e-8),1e-8)),np.log(np.maximum(delta,1e-8))]
        sx=px[:,:,None]+px[:,None,:]; sy=py[:,:,None]+py[:,None,:]; sz=pz[:,:,None]+pz[:,None,:]; se=e[:,:,None]+e[:,None,:]
        f.append(np.log(np.maximum(se*se-sx*sx-sy*sy-sz*sz,1e-8)))
        z=np.stack(f,-1)
        z=self._batchnorm(z,"pair_embed.embed.0")
        for conv,bni,last in [("pair_embed.embed.1","pair_embed.embed.2",False),("pair_embed.embed.4","pair_embed.embed.5",False),("pair_embed.embed.7","pair_embed.embed.8",False),("pair_embed.embed.10","pair_embed.embed.11",True)]:
            z=self._linear(z,conv); z=self._batchnorm(z,bni)
            if not last or self.pair_final_gelu:z=gelu(z)
        return z.transpose(0,3,1,2)
    def block(self,x,prefix,bias=None,cls=False,mask=None):
        if cls:
            u=np.concatenate([x[0],x[1]],axis=1); query=x[0]; residual=query; z=self._norm(u,prefix+".pre_attn_norm")
            qkv_w=self.s[prefix+".attn.in_proj.weight"];qkv_b=self.s[prefix+".attn.in_proj.bias"]
            q=self.lin(query,qkv_w[:128],qkv_b[:128]); k=self.lin(z,qkv_w[128:256],qkv_b[128:256]);v=self.lin(z,qkv_w[256:],qkv_b[256:])
        else:
            residual=x
            z=self._norm(x,prefix+".pre_attn_norm"); qkv=self.lin(z,self.s[prefix+".attn.in_proj.weight"],self.s[prefix+".attn.in_proj.bias"]);q,k,v=np.split(qkv,3,-1)
        B,Nq,_=q.shape;Nk=k.shape[1];H=8;D=16
        q=q.reshape(B,Nq,H,D).transpose(0,2,1,3);k=k.reshape(B,Nk,H,D).transpose(0,2,1,3);v=v.reshape(B,Nk,H,D).transpose(0,2,1,3)
        score=q@k.transpose(0,1,3,2)/4.0
        if bias is not None:score+=bias
        if mask is not None: score=np.where(mask[:,None,None,:],score,-np.inf)
        z=(softmax(score)@v).transpose(0,2,1,3).reshape(B,Nq,128)
        z=self._linear(z,prefix+".attn.out_proj").reshape(B,Nq,H,D)
        z*=self.s[prefix+".c_attn"].reshape(1,1,H,1)
        if self.model_variant=="legacy_mha":z=z.transpose(0,1,3,2)
        z=self._norm(z.reshape(B,Nq,128),prefix+".post_attn_norm");x=residual+z
        z=self._norm(x,prefix+".pre_fc_norm");z=gelu(self._linear(z,prefix+".fc1"));z=self._norm(z,prefix+".post_fc_norm");z=self._linear(z,prefix+".fc2")
        return self.s[prefix+".w_resid"]*x+z
    def __call__(self,x,v,mask):
        # x B,C,N -> token-major; mask B,N bool
        z=x.transpose(0,2,1);z=self._batchnorm(z,"embed.input_bn")
        for i in (0,3,6):
            z=self._norm(z,f"embed.embed.{i}");z=self._linear(z,f"embed.embed.{i+1}")
            z=gelu(z)
        z=np.where(mask[:,:,None],z,0)
        pair=self.pair(v)
        pair=np.where((mask[:,None,:,None]&mask[:,None,None,:]),pair,0)
        for i in range(8):z=self.block(z,f"blocks.{i}",pair,mask=mask)
        cls=np.broadcast_to(self.s["cls_token"],(z.shape[0],1,128)).copy(); cmask=np.concatenate([np.ones((z.shape[0],1),bool),mask],1)
        for i in range(2):cls=self.block((cls,z),f"cls_blocks.{i}",cls=True,mask=cmask)
        z=self._norm(cls[:,0],"norm");return self._linear(z,"fc.0")
