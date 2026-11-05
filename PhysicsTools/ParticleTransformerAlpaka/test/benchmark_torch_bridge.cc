// C++ TorchScript forward bridge: same LibTorch inference engine as PyTorchAlpaka,
// without Python dispatch in the timed batch loop. CMSSW wrapping is excluded.
#include <torch/script.h>
#include <c10/core/InferenceMode.h>
#include <ATen/Context.h>
#include <string>
#include <memory>

struct TorchHandle {
  torch::jit::Module model;
  torch::Tensor features,vectors,mask,points,output;
  int jets;
  TorchHandle(torch::jit::Module m):model(std::move(m)){}
};
static thread_local std::string lastError;
extern "C" const char* ts_error(){return lastError.c_str();}
extern "C" void* ts_create(const char* path,int device){
 try {
  c10::InferenceMode guard;
  auto h=std::make_unique<TorchHandle>(torch::jit::load(path,torch::Device(torch::kCUDA,device)));
  h->model.eval();return h.release();
 }catch(std::exception const&e){lastError=e.what();return nullptr;}
}
extern "C" int ts_bind(void* p,float* f,float* v,float* m,float* points,float* output,int jets,int device,int soa){
 try {
  c10::InferenceMode guard;auto& h=*static_cast<TorchHandle*>(p);
  auto opts=torch::TensorOptions().device(torch::Device(torch::kCUDA,device)).dtype(torch::kFloat32);
  int64_t pitch=((jets+31)/32)*32;
  auto stride3=std::vector<int64_t>{1,pitch,16*pitch};
  auto stride2=std::vector<int64_t>{1,pitch};
  h.features=soa ? torch::from_blob(f,{jets,16,15},stride3,opts) : torch::from_blob(f,{jets,16,15},opts);
  h.vectors=soa ? torch::from_blob(v,{jets,16,4},stride3,opts) : torch::from_blob(v,{jets,16,4},opts);
  h.mask=soa ? torch::from_blob(m,{jets,16},stride2,opts) : torch::from_blob(m,{jets,16},opts);
  h.points=soa ? torch::from_blob(points,{jets,16,2},stride3,opts) : torch::from_blob(points,{jets,16,2},opts);
  h.output=torch::from_blob(output,{jets,1},opts);h.jets=jets;return 0;
 }catch(std::exception const&e){lastError=e.what();return -1;}
}
extern "C" int ts_run(void* p,int batch){
 try {
  c10::InferenceMode guard;auto& h=*static_cast<TorchHandle*>(p);
  for(int first=0;first<h.jets;first+=batch){
   int end=std::min(first+batch,h.jets);
   auto result=h.model.forward({h.points.slice(0,first,end),h.features.slice(0,first,end),
      h.vectors.slice(0,first,end),h.mask.slice(0,first,end)}).toTensor();
   h.output.slice(0,first,end).copy_(result.reshape({end-first,1}));
  }return 0;
 }catch(std::exception const&e){lastError=e.what();return -1;}
}
extern "C" void ts_destroy(void* p){delete static_cast<TorchHandle*>(p);}
