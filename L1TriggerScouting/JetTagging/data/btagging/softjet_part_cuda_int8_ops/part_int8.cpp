
#include <ATen/ATen.h>
#include <torch/library.h>
#include <optional>
#include <vector>

void check_inputs(const at::Tensor& x, const at::Tensor& w,
                  const at::Tensor& s, const std::optional<at::Tensor>& b) {
  TORCH_CHECK(x.scalar_type() == at::kFloat, "part_int8: input must be float32");
  TORCH_CHECK(w.scalar_type() == at::kChar, "part_int8: weights must be int8");
  TORCH_CHECK(s.scalar_type() == at::kFloat, "part_int8: scales must be float32");
  TORCH_CHECK(x.layout() == at::kStrided && w.layout() == at::kStrided &&
              s.layout() == at::kStrided, "part_int8: expected strided tensors");
  TORCH_CHECK(x.dim() >= 1 && w.dim() == 2 && s.dim() == 1,
              "part_int8: expected input [..., K], weights [K, N], scales [N]");
  TORCH_CHECK(x.size(-1) == w.size(0) && s.size(0) == w.size(1),
              "part_int8: inconsistent tensor dimensions");
  TORCH_CHECK(w.size(0) > 0 && w.size(1) > 0, "part_int8: empty weight matrix");
  TORCH_CHECK(x.device() == w.device() && x.device() == s.device(),
              "part_int8: tensors must be on the same device");
  if (b.has_value()) {
    TORCH_CHECK(b->scalar_type() == at::kFloat && b->dim() == 1 &&
                b->size(0) == w.size(1) && b->device() == x.device() &&
                b->layout() == at::kStrided, "part_int8: invalid bias");
  }
}

at::Tensor linear_cpu(const at::Tensor& x, const at::Tensor& w,
                      const at::Tensor& s, const std::optional<at::Tensor>& b) {
  check_inputs(x, w, s, b);
  auto shape = x.sizes().vec();
  shape.back() = w.size(1);
  auto result = at::matmul(x.reshape({-1, w.size(0)}), w.to(at::kFloat) * s);
  if (b.has_value()) result.add_(*b);
  return result.reshape(shape);
}

#ifdef PART_INT8_WITH_CUDA
at::Tensor linear_cuda(const at::Tensor&, const at::Tensor&, const at::Tensor&,
                        const std::optional<at::Tensor>&);
#endif

TORCH_LIBRARY(part_int8, m) {
  m.def("linear(Tensor input, Tensor weight, Tensor scale, Tensor? bias=None) -> Tensor");
}
TORCH_LIBRARY_IMPL(part_int8, CPU, m) { m.impl("linear", &linear_cpu); }
#ifdef PART_INT8_WITH_CUDA
TORCH_LIBRARY_IMPL(part_int8, CUDA, m) { m.impl("linear", &linear_cuda); }
#endif
