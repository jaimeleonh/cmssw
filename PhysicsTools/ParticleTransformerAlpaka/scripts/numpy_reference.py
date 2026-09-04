"""NumPy reconstruction of the exported graph, built only from the emitted
header and binary blobs.

This is the independent check behind VALIDATION.md: it shares no code with the
C++ kernels, so agreement between the two means the generated descriptors, the
tensor ordering and the offsets are all consistent.  Point it at the package's
own data directory and it reproduces the scalar reference kernel to 1e-9.
"""
import numpy as np, re
from scipy.special import erf
import pathlib
_pkg=pathlib.Path(__file__).parents[1]
hdr=(_pkg/'interface/alpaka/GeneratedModel.h').read_text()
W=np.fromfile(_pkg/'data/model_weights.int8.bin',np.int8); P=np.fromfile(_pkg/'data/model_params.fp32.bin',np.float32)
L=re.findall(r'Linear\{(\d+),(\d+),m\.weights\+(\d+),m\.params\+(\d+),m\.params\+(\d+),m\.params\+(\d+)\}',hdr)
N=re.findall(r'Norm\{(\d+),1e-5f,m\.params\+(\d+),m\.params\+(\d+)\}',hdr)
BN=re.findall(r'BatchNorm\{(\d+),1e-5f,m\.params\+(\d+),m\.params\+(\d+),m\.params\+(\d+),m\.params\+(\d+)\}',hdr)
def blockscales(fn,cnt):
    out=[]
    body=re.search(rf'inline BlockData {fn}\(ModelView m,int i\)\{{switch\(i\)\{{(.*?)default',hdr,re.S).group(1)
    for part in re.findall(r'BlockData\{(.*?)\}\;',body+';'):
        o=[int(v) for v in re.findall(r'm\.params\+(\d+)',part)[-2:]]; out.append(o)
    return out[:cnt]
CA=blockscales('block',8); CAC=blockscales('clsBlock',2)
cls_tok=int(re.search(r'clsToken\(ModelView m\)\{return m\.params\+(\d+)',hdr).group(1))
def lin(i,x):
    K,O,w,sc,b,a=[int(v) for v in L[i]]
    q=W[w:w+K*O].reshape(O,K).astype(np.float64)*P[sc:sc+O][:,None].astype(np.float64)
    return x@q.T+P[b:b+O]
def nrm(i,x):
    S,g,be=[int(v) for v in N[i]]
    return (x-x.mean(-1,keepdims=True))/np.sqrt(x.var(-1,keepdims=True)+1e-5)*P[g:g+S]+P[be:be+S]
def bn(i,x):
    S,g,be,rm,rv=[int(v) for v in BN[i]]
    return (x-P[rm:rm+S])/np.sqrt(P[rv:rv+S]+1e-5)*P[g:g+S]+P[be:be+S]
gelu=lambda x:.5*x*(1+erf(x/np.sqrt(2)))

def attend(Q,Kk,V,valid,bias=None):
    T,_=Q.shape; out=np.zeros((Q.shape[0],128))
    for h in range(8):
        q=Q[:,h*16:(h+1)*16]; k=Kk[:,h*16:(h+1)*16]; v=V[:,h*16:(h+1)*16]
        s=0.25*(q@k.T)
        if bias is not None: s=s+bias[h]
        s=np.where(valid[None,:],s,-3.4e38); s=s-s.max(-1,keepdims=True)
        e=np.exp(s); out[:,h*16:(h+1)*16]=(e/e.sum(-1,keepdims=True))@v
    return out

def run(feat,mask,pair=None):
    h=bn(0,feat.astype(np.float64))
    h=gelu(lin(0,nrm(0,h))); h=gelu(lin(1,nrm(1,h))); h=gelu(lin(2,nrm(2,h)))
    x=np.where(mask[:,None],h,0)
    for b in range(8):
        ni,li=3+4*b,7+4*b
        qkv=lin(li,nrm(ni,x)); Q,Kk,V=qkv[:,:128],qkv[:,128:256],qkv[:,256:]
        o=lin(li+1,attend(Q,Kk,V,mask,pair))
        o=(o.reshape(-1,8,16)*P[CA[b][0]:CA[b][0]+8][None,:,None]).reshape(-1,128)
        x=x+nrm(ni+1,o); r=x
        f=gelu(lin(li+2,nrm(ni+2,x))); f=lin(li+3,nrm(ni+3,f))
        x=P[CA[b][1]:CA[b][1]+128]*r+f
    c=P[cls_tok:cls_tok+128].astype(np.float64)[None,:]
    for b in range(2):
        ni,li=3+32+4*b,7+32+4*b
        u=nrm(ni,np.concatenate([c,x],0))
        qkv_q=lin(li,c); Q=qkv_q[:,:128]
        qkv_u=lin(li,u); Kk,V=qkv_u[:,128:256],qkv_u[:,256:]
        valid=np.concatenate([[True],mask])
        o=lin(li+1,attend(Q,Kk,V,valid))
        o=(o.reshape(-1,8,16)*P[CAC[b][0]:CAC[b][0]+8][None,:,None]).reshape(-1,128)
        c=c+nrm(ni+1,o); r=c
        f=gelu(lin(li+2,nrm(ni+2,c))); f=lin(li+3,nrm(ni+3,f))
        c=P[CAC[b][1]:CAC[b][1]+128]*r+f
    z=nrm(3+40,c)
    fw,fb=[int(v) for v in re.search(r'FloatLinear\{128,2,m\.params\+(\d+),m\.params\+(\d+)\}',hdr).groups()]
    lg=(z@P[fw:fw+256].reshape(2,128).T.astype(np.float64)+P[fb:fb+2])[0]
    e=np.exp(lg-lg.max()); return e/e.sum()
