// Generated from an INT8 QAT DeepSets checkpoint (deepset_int8.py); values live in external binary files; do not edit.
#pragma once
#include <cstddef>
#include <cstdint>
namespace generated {
inline constexpr char checkpoint_sha256[]="eb1884c55a2a9f3fd0e4a1cb9dfc3bd28e113d352d68a31f2f3ef54c08752d77";
inline constexpr char architecture[]="deepsets_int8_v1";
inline constexpr char storage_precision[]="int8_weights_fp32_params_with_activation_scales";
inline constexpr char weights_sha256[]="23701afaf4ad5499cdf995f3c6f93e7d719b14cfe0cdd760491cb353667dc2fd";
inline constexpr char params_sha256[]="650cc09c22c177569945e19e4d27df15c6ab3b707f7f204f4e276a3128f47fa5";
inline constexpr std::size_t weights_size=2048;
inline constexpr std::size_t params_size=1096;
inline constexpr int input_features=15;
inline constexpr int embed_width=32;
inline constexpr int num_classes=2;
struct ModelView { int8_t const* weights; float const* params; };
ALPAKA_FN_ACC inline FloatLinear inputLinear(ModelView m){return FloatLinear{15,32,m.params+0,m.params+480};}
ALPAKA_FN_ACC inline BatchNorm phiNorm(ModelView m,int i){switch(i){case 0:return BatchNorm{32,1e-5f,m.params+512,m.params+544,m.params+576,m.params+608};case 1:return BatchNorm{32,1e-5f,m.params+707,m.params+739,m.params+771,m.params+803};default:return {};}}
ALPAKA_FN_ACC inline Linear phiProject(ModelView m){return Linear{32,32,m.weights+0,m.params+640,m.params+672,m.params+704};}
ALPAKA_FN_ACC inline Linear rhoProject(ModelView m){return Linear{32,32,m.weights+1024,m.params+835,m.params+867,m.params+899};}
ALPAKA_FN_ACC inline BatchNorm rhoNorm(ModelView m){return BatchNorm{32,1e-5f,m.params+902,m.params+934,m.params+966,m.params+998};}
ALPAKA_FN_ACC inline FloatLinear classifier(ModelView m){return FloatLinear{32,2,m.params+1030,m.params+1094};}
} // namespace generated
