import hashlib
import pathlib

import numpy as np

def test_per_channel_int8_error_bound():
    rng=np.random.default_rng(7); w=rng.normal(size=(37,53)).astype("f4")
    s=np.maximum(np.max(np.abs(w),axis=1),1e-12)/127
    q=np.clip(np.rint(w/s[:,None]),-127,127).astype("i1")
    err=np.abs(w-q.astype("f4")*s[:,None])
    assert np.all(err <= s[:,None]/2 + 2e-7)


def test_external_model_files():
    package = pathlib.Path(__file__).parents[1]
    header = (package / "interface/alpaka/GeneratedModel.h").read_text()
    weights = (package / "data/model_weights.int8.bin").read_bytes()
    params = (package / "data/model_params.fp32.bin").read_bytes()
    # 2,108,032 INT8 weights: every dense tensor except the classifier, which
    # quantization-aware training kept in FP32 and which therefore lives in the
    # parameter blob instead.
    assert len(weights) == 2_108_032
    # 48 layers x 3 activation scales on top of the previous parameter blob.
    assert len(params) == 59_713 * 4
    assert "FloatLinear classifier" in header
    assert "host_weights" not in header and "host_params" not in header
    assert hashlib.sha256(weights).hexdigest() in header
    assert hashlib.sha256(params).hexdigest() in header
