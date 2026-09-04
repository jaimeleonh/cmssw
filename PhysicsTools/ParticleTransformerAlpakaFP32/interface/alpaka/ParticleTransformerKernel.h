#ifndef PhysicsTools_ParticleTransformerAlpakaFP32_ParticleTransformerKernel_h
#define PhysicsTools_ParticleTransformerAlpakaFP32_ParticleTransformerKernel_h

// Scalar, single-threaded reference evaluation of the exported graph.
//
// This is not the production path: CMSSW runs the block-cooperative kernel in
// ParticleTransformerFastKernel.h.  The reference is kept because it is the
// simplest possible transcription of the Weaver graph and is what the
// regression test compares the fast kernel against.  It reads the weights in
// their original row-major [out][in] layout.

#include <alpaka/alpaka.hpp>
#include <cstdint>

#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ModelDefinition.h"

namespace partfp32 {
template <typename A> ALPAKA_FN_ACC inline void linear(A const&,generated::Linear const& l,float const* x,float* y){for(int o=0;o<l.out;++o){float z=l.bias[o];for(int i=0;i<l.in;++i)z+=x[i]*l.weight[o*l.in+i];y[o]=z;}}
template <typename A> ALPAKA_FN_ACC inline void norm(A const& acc,generated::Norm const& p,float const* x,float* y){float m=0,v=0;for(int i=0;i<p.size;++i)m+=x[i];m/=p.size;for(int i=0;i<p.size;++i){float d=x[i]-m;v+=d*d;}v/=p.size;float r=1.f/alpaka::math::sqrt(acc,v+p.eps);for(int i=0;i<p.size;++i)y[i]=(x[i]-m)*r*p.gamma[i]+p.beta[i];}
template <typename A> ALPAKA_FN_ACC inline void batchNorm(A const& acc,generated::BatchNorm const& p,float const* x,float* y){for(int i=0;i<p.size;++i)y[i]=(x[i]-p.mean[i])/alpaka::math::sqrt(acc,p.var[i]+p.eps)*p.gamma[i]+p.beta[i];}
template <typename A> ALPAKA_FN_ACC inline float gelu(A const& acc,float x){return .5f*x*(1.f+alpaka::math::erf(acc,x*.7071067811865475f));}

template <typename A> ALPAKA_FN_ACC inline void encodePairs(A const& acc,generated::ModelView model,float const* p4,bool const* mask,float* pair){
  constexpr int N=16;float z0[4],z1[64],z2[64];
  for(int i=0;i<N;++i)for(int j=0;j<N;++j){float pxi=p4[i*4],pyi=p4[i*4+1],pzi=p4[i*4+2],ei=p4[i*4+3],pxj=p4[j*4],pyj=p4[j*4+1],pzj=p4[j*4+2],ej=p4[j*4+3];float pti=alpaka::math::sqrt(acc,pxi*pxi+pyi*pyi),ptj=alpaka::math::sqrt(acc,pxj*pxj+pyj*pyj);float ri=.5f*alpaka::math::log(acc,1.f+2.f*pzi/alpaka::math::max(acc,ei-pzi,1e-20f)),rj=.5f*alpaka::math::log(acc,1.f+2.f*pzj/alpaka::math::max(acc,ej-pzj,1e-20f));float pi=alpaka::math::atan2(acc,pyi,pxi),pj=alpaka::math::atan2(acc,pyj,pxj),dp=pi-pj;constexpr float twoPi=6.283185307179586f;while(dp>3.141592653589793f)dp-=twoPi;while(dp<=-3.141592653589793f)dp+=twoPi;float d=alpaka::math::sqrt(acc,(ri-rj)*(ri-rj)+dp*dp),ptm=alpaka::math::min(acc,pti,ptj);z0[0]=alpaka::math::log(acc,alpaka::math::max(acc,ptm*d,1e-8f));z0[1]=alpaka::math::log(acc,alpaka::math::max(acc,ptm/alpaka::math::max(acc,pti+ptj,1e-8f),1e-8f));z0[2]=alpaka::math::log(acc,alpaka::math::max(acc,d,1e-8f));float sx=pxi+pxj,sy=pyi+pyj,sz=pzi+pzj,se=ei+ej;z0[3]=alpaka::math::log(acc,alpaka::math::max(acc,se*se-sx*sx-sy*sy-sz*sz,1e-8f));auto pbn=generated::pairInputBN(model);batchNorm(acc,pbn,z0,z1);float* cur=z1;float* nxt=z2;for(int l=0;l<4;++l){auto pl=generated::pairLinear(model,l);auto bn=generated::pairBN(model,l);linear(acc,pl,cur,nxt);batchNorm(acc,bn,nxt,cur);if(l<3||generated::pair_final_gelu)for(int q=0;q<bn.size;++q)cur[q]=gelu(acc,cur[q]);}for(int h=0;h<8;++h)pair[h*N*N+i*N+j]=(mask[i]&&mask[j])?cur[h]:0.f;}
}

template <typename A> ALPAKA_FN_ACC inline void transformerBlock(A const& acc,generated::BlockData const& p,float* x,bool const* mask,float const* pair){
  constexpr int N=16,D=128,H=8,HD=16;float nrm[N*D],qkv[N*3*D],attn[N*D],hidden[N*512],work[512],score[N];for(int t=0;t<N;++t){norm(acc,p.preAttn,&x[t*D],&nrm[t*D]);linear(acc,p.qkv,&nrm[t*D],&qkv[t*3*D]);}for(int tq=0;tq<N;++tq){for(int h=0;h<H;++h){float mx=-3.4e38f;for(int tk=0;tk<N;++tk){float s=0;for(int d=0;d<HD;++d)s+=qkv[tq*3*D+h*HD+d]*qkv[tk*3*D+D+h*HD+d];s=.25f*s+pair[h*N*N+tq*N+tk];if(!mask[tk])s=-3.4e38f;score[tk]=s;if(s>mx)mx=s;}float den=0;for(int tk=0;tk<N;++tk){score[tk]=alpaka::math::exp(acc,score[tk]-mx);den+=score[tk];}for(int d=0;d<HD;++d){float s=0;for(int tk=0;tk<N;++tk)s+=score[tk]/den*qkv[tk*3*D+2*D+h*HD+d];attn[tq*D+h*HD+d]=s;}}linear(acc,p.out,&attn[tq*D],work);for(int h=0;h<H;++h)for(int d=0;d<HD;++d){int dst=generated::transpose_scaled_heads?d*H+h:h*HD+d;attn[tq*D+dst]=work[h*HD+d]*p.headScale[h];}norm(acc,p.postAttn,&attn[tq*D],work);for(int d=0;d<D;++d)x[tq*D+d]+=work[d];}for(int t=0;t<N;++t){norm(acc,p.preFc,&x[t*D],work);linear(acc,p.fc1,work,&hidden[t*512]);for(int d=0;d<512;++d)hidden[t*512+d]=gelu(acc,hidden[t*512+d]);norm(acc,p.postFc,&hidden[t*512],work);linear(acc,p.fc2,work,&attn[t*D]);for(int d=0;d<D;++d)x[t*D+d]=p.residualScale[d]*x[t*D+d]+attn[t*D+d];}
}

template <typename A> ALPAKA_FN_ACC inline void classBlock(A const& acc,generated::BlockData const& p,float* cls,float const* x,bool const* mask){
  constexpr int N=16,D=128,H=8,HD=16;float u[17*D],nu[17*D],qkv[17*3*D],query[D],merged[D],work[512],score[17];for(int d=0;d<D;++d)u[d]=cls[d];for(int t=0;t<N;++t)for(int d=0;d<D;++d)u[(t+1)*D+d]=x[t*D+d];for(int t=0;t<17;++t){norm(acc,p.preAttn,&u[t*D],&nu[t*D]);linear(acc,p.qkv,&nu[t*D],&qkv[t*3*D]);}for(int o=0;o<D;++o){float value=p.qkv.bias[o];for(int i=0;i<D;++i)value+=cls[i]*p.qkv.weight[o*D+i];query[o]=value;}for(int h=0;h<H;++h){float mx=-3.4e38f;for(int tk=0;tk<17;++tk){float s=0;for(int d=0;d<HD;++d)s+=query[h*HD+d]*qkv[tk*3*D+D+h*HD+d];s*=.25f;if(tk>0&&!mask[tk-1])s=-3.4e38f;score[tk]=s;if(s>mx)mx=s;}float den=0;for(int tk=0;tk<17;++tk){score[tk]=alpaka::math::exp(acc,score[tk]-mx);den+=score[tk];}for(int d=0;d<HD;++d){float s=0;for(int tk=0;tk<17;++tk)s+=score[tk]/den*qkv[tk*3*D+2*D+h*HD+d];merged[h*HD+d]=s;}}linear(acc,p.out,merged,work);for(int h=0;h<H;++h)for(int d=0;d<HD;++d){int dst=generated::transpose_scaled_heads?d*H+h:h*HD+d;merged[dst]=work[h*HD+d]*p.headScale[h];}norm(acc,p.postAttn,merged,work);for(int d=0;d<D;++d)cls[d]+=work[d];norm(acc,p.preFc,cls,merged);linear(acc,p.fc1,merged,work);for(int d=0;d<512;++d)work[d]=gelu(acc,work[d]);norm(acc,p.postFc,work,u);linear(acc,p.fc2,u,merged);for(int d=0;d<D;++d)cls[d]=p.residualScale[d]*cls[d]+merged[d];
}

struct Kernel{template <typename A> ALPAKA_FN_ACC void operator()(A const& acc,generated::ModelView model,float const* features,float const* p4,bool const* mask,float* output)const{constexpr int N=16,D=128;float x[N*D],a[512],b[512],pair[8*N*N],cls[D];auto ibn=generated::inputBN(model);for(int t=0;t<N;++t){batchNorm(acc,ibn,&features[t*15],a);for(int l=0;l<3;++l){auto en=generated::embedNorm(model,l);auto el=generated::embedLinear(model,l);norm(acc,en,a,b);linear(acc,el,b,a);for(int d=0;d<el.out;++d)a[d]=gelu(acc,a[d]);}for(int d=0;d<D;++d)x[t*D+d]=mask[t]?a[d]:0.f;}encodePairs(acc,model,p4,mask,pair);for(int l=0;l<8;++l)transformerBlock(acc,generated::block(model,l),x,mask,pair);auto ct=generated::clsToken(model);for(int d=0;d<D;++d)cls[d]=ct[d];for(int l=0;l<2;++l)classBlock(acc,generated::clsBlock(model,l),cls,x,mask);norm(acc,generated::finalNorm(model),cls,a);linear(acc,generated::classifier(model),a,b);float mx=b[0]>b[1]?b[0]:b[1],e0=alpaka::math::exp(acc,b[0]-mx),e1=alpaka::math::exp(acc,b[1]-mx);output[0]=e0/(e0+e1);output[1]=e1/(e0+e1);}};
}
#endif
