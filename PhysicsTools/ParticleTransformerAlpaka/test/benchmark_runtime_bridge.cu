// Standalone CUDA bridge to the production graph interpreter.
// Used by benchmark_runtime_comparison.py; excludes CMSSW event/SoA overhead.
#include <cuda_runtime.h>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/ModelLoader.h"

struct Acc {};
struct GraphHandle {
  nn::ModelView view{};
  nn::OpRecord* ops = nullptr;
  nn::BufferRecord* buffers = nullptr;
  std::int8_t* weights = nullptr;
  float* params = nullptr;
  bool fitsWarp = false;
  std::size_t bytes = 0;
};
static thread_local std::string lastError;
static void check(cudaError_t code) {
  if (code != cudaSuccess) throw std::runtime_error(cudaGetErrorString(code));
}
template<class T> static T* upload(T const* src, std::size_t count) {
  if (!count) return nullptr;
  T* dst = nullptr;
  check(cudaMalloc(&dst, count * sizeof(T)));
  check(cudaMemcpy(dst, src, count * sizeof(T), cudaMemcpyHostToDevice));
  return dst;
}
template<class State> __device__ void stage(nn::ModelView const& model, State& s,
    float const* features, float const* vectors, float const* masks, int jet, int lane, int lanes, int jets, bool soa) {
  int n = model.maxTokens, f = model.features.cols;
  int pitch = ((jets+31)/32)*32;
  for (int at = lane; at < n*f; at += lanes)
    s.arena[model.features.offset+(at/f)*model.features.stride+at%f] = features[soa ? (long)jet+(at/f)*(long)pitch+(at%f)*(long)n*pitch : (long)jet*n*f+at];
  if (model.vectors.valid()) {
    int c = model.vectors.cols;
    for (int at = lane; at < n*c; at += lanes)
      s.arena[model.vectors.offset+(at/c)*model.vectors.stride+at%c] = vectors[soa ? (long)jet+(at/c)*(long)pitch+(at%c)*(long)n*pitch : (long)jet*n*c+at];
  }
  for (int at = lane; at < n; at += lanes) s.mask[at] = masks[soa ? (long)jet+(long)at*pitch : (long)jet*n+at];
}
__global__ __launch_bounds__(nn::kThreadsPerJet) void graphBlock(nn::ModelView model,
    float const* features, float const* vectors, float const* masks, float* output, int jets, bool soa) {
  __shared__ nn::GraphShared shared;
  Acc acc;
  for (int jet = blockIdx.x; jet < jets; jet += gridDim.x) {
    stage(model, shared, features, vectors, masks, jet, threadIdx.x, blockDim.x, jets, soa);
    __syncthreads();
    nn::runGraph(acc, model, shared, nn::GraphInputs{nullptr,nullptr});
    if (threadIdx.x == 0) output[jet] = shared.output[0];
    __syncthreads();
  }
}
__global__ __launch_bounds__(nn::kThreadsPerJet) void graphWarp(nn::ModelView model,
    float const* features, float const* vectors, float const* masks, float* output, int jets, bool soa) {
  __shared__ nn::GraphWarpBlock states;
  Acc acc;
  int warp=threadIdx.x/nn::kWarpLanes, lane=threadIdx.x%nn::kWarpLanes;
  nn::SubGroup<Acc> group{acc,lane,nn::kWarpLanes};
  auto& shared=states.jet[warp];
  for(int jet=blockIdx.x*nn::kWarpJetsPerBlock+warp;jet<jets;jet+=gridDim.x*nn::kWarpJetsPerBlock) {
    stage(model,shared,features,vectors,masks,jet,lane,nn::kWarpLanes,jets,soa);
    nn::barrier(group);
    nn::runGraph(group,model,shared,nn::GraphInputs{nullptr,nullptr});
    if(lane==0) output[jet]=shared.output[0];
    nn::barrier(group);
  }
}
extern "C" const char* graph_error() {return lastError.c_str();}
extern "C" void graph_destroy(void* pointer) {
  auto* h=static_cast<GraphHandle*>(pointer); if(!h)return;
  cudaFree(h->ops);cudaFree(h->buffers);cudaFree(h->weights);cudaFree(h->params);delete h;
}
extern "C" void* graph_create(const char* path) {
  GraphHandle* h=nullptr;
  try {
    std::ifstream file(path,std::ios::binary|std::ios::ate);
    if(!file)throw std::runtime_error("Cannot open graph");
    std::vector<char> image(file.tellg());file.seekg(0);file.read(image.data(),image.size());
    if(!file)throw std::runtime_error("Cannot read graph");
    nn::LoadedModel m;std::string error;
    if(!nn::loadModel(image.data(),image.size(),m,error))throw std::runtime_error(error);
    h=new GraphHandle;
    h->ops=upload(m.ops,m.header->opCount);h->buffers=upload(m.buffers,m.header->bufferCount);
    h->weights=upload(m.weights,m.header->weightsCount);h->params=upload(m.params,m.header->paramsCount);
    h->view=m.view(h->weights,h->params,h->ops,h->buffers);
    h->fitsWarp=m.fitsWarpMode();
    h->bytes=m.header->opCount*sizeof(nn::OpRecord)+m.header->bufferCount*sizeof(nn::BufferRecord)
      +m.header->weightsCount+m.header->paramsCount*sizeof(float);
    return h;
  } catch(std::exception const& e){lastError=e.what();graph_destroy(h);return nullptr;}
}
extern "C" std::size_t graph_bytes(void* p){return static_cast<GraphHandle*>(p)->bytes;}
extern "C" int graph_fits_warp(void* p){return static_cast<GraphHandle*>(p)->fitsWarp;}
extern "C" int graph_launch(void* p,float const* f,float const* v,float const* m,float* o,
    int jets,int warp,int soa,void* stream) {
  auto* h=static_cast<GraphHandle*>(p);
  if(warp&&!h->fitsWarp){lastError="Graph does not fit warp arena";return -1;}
  auto s=reinterpret_cast<cudaStream_t>(stream);
  if(warp) graphWarp<<<(jets+nn::kWarpJetsPerBlock-1)/nn::kWarpJetsPerBlock,nn::kThreadsPerJet,0,s>>>(h->view,f,v,m,o,jets,soa);
  else graphBlock<<<jets,nn::kThreadsPerJet,0,s>>>(h->view,f,v,m,o,jets,soa);
  auto e=cudaGetLastError();if(e!=cudaSuccess){lastError=cudaGetErrorString(e);return -1;}return 0;
}
extern "C" int graph_resources(int warp, long long* result) {
  try {
    cudaFuncAttributes a{};int blocks=0;
    if(warp){check(cudaFuncGetAttributes(&a,graphWarp));check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&blocks,graphWarp,nn::kThreadsPerJet,0));}
    else{check(cudaFuncGetAttributes(&a,graphBlock));check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&blocks,graphBlock,nn::kThreadsPerJet,0));}
    result[0]=a.sharedSizeBytes;result[1]=a.localSizeBytes;result[2]=a.numRegs;result[3]=blocks;
    return 0;
  }catch(std::exception const&e){lastError=e.what();return -1;}
}
