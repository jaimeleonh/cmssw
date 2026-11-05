#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "scripts"))
from weaver_checkpoint import detect_checkpoint_variant, load_state_dict
from weaver_reference import WeaverParT, softmax


PRINTED_JET = np.array(
    [
        [0.013672,-0.013212,1.202409,1.105805,2.427650,2.430056,-0.723950,0.013672,-0.013212,0.104687,0.062500,1,1,0,0,0,0,-28.023132,-12.039688,18.866461,35.863541],
        [-0.008145,0.026058,0.728588,0.624072,1.953830,1.948323,-0.690794,-0.008145,0.026058,-0.045313,0.062500,1,1,0,0,0,0,-13.990071,-6.672922,9.192499,18.020878],
        [0.052942,-0.021938,0.717206,0.635448,1.942447,1.959699,-0.570771,0.052942,-0.021938,0.054688,0.062500,1,1,0,0,0,0,-14.063565,-5.897343,10.144864,18.316133],
        [-0.112864,-0.026302,0.549435,0.410465,1.774676,1.734716,-0.336446,-0.112864,-0.026302,-0.045313,0.093750,1,1,0,0,0,0,-11.086554,-4.592202,5.692147,13.281587],
        [0.039852,0.056602,0.220432,0.133641,1.445674,1.457892,-0.523106,0.039852,0.056602,0.054688,0.109375,1,1,0,0,0,0,-6.667628,-3.434055,4.871789,8.943395],
        [-0.095411,-0.109205,0.092807,-0.040840,1.318049,1.283411,-0.219945,-0.095411,-0.109205,0.204687,0.203125,1,1,0,0,0,0,-5.952474,-1.905402,3.085850,6.970292],
        [0.096575,0.165684,-0.219594,-0.283980,1.005648,1.040271,-0.032896,0.096575,0.165684,-0.145313,0.125000,1,1,0,0,0,0,-3.335543,-2.207748,2.873172,4.924949],
        [-0.021235,0.122051,-0.481879,-0.591025,0.743362,0.733226,-0.304461,-0.021235,0.122051,-0.145313,0.187500,1,1,0,0,0,0,-2.357210,-1.416355,1.589214,3.176177],
        [-0.191404,0.100235,-0.548597,-0.709327,0.676645,0.614924,0.064246,-0.191404,0.100235,-0.145313,0.203125,1,1,0,0,0,0,-2.170497,-1.240541,0.971980,2.682302],
        [-0.169588,-0.017575,-0.622349,-0.777402,0.602893,0.546849,-0.118017,-0.169588,-0.017575,0.054688,0.203125,1,1,0,0,0,0,-2.071136,-0.879145,0.927662,2.433733],
        [0.083485,0.095871,-0.622349,-0.692041,0.602893,0.632210,-0.291495,0.083485,0.095871,2.304688,0,0.949219,0,0,1,0,0,-1.958300,-1.107953,1.580034,2.749365],
        [-0.012508,-0.187745,-0.906174,-1.012244,0.319067,0.312007,-0.047357,-0.012508,-0.187745,2.304688,0,0.292969,0,1,0,0,0,-1.460069,-0.343801,0.881996,1.740091],
        [0.275471,-0.109205,-1.033800,-1.018251,0.191442,0.306001,0.385311,0.275471,-0.109205,2.304688,0,0.183594,0,1,0,0,0,-1.190495,-0.381080,1.189073,1.725223],
        [0.240565,-0.078661,-1.033800,-1.034864,0.191442,0.289387,0.212395,0.240565,-0.078661,2.304688,0,0.839844,0,0,1,0,0,-1.178302,-0.417259,1.129563,1.684759],
    ], dtype=np.float32)


def main() -> None:
    package = pathlib.Path(__file__).parents[1]
    root = package.parent
    state = load_state_dict(root / "upload/ParT_checkpoint.pt")
    variant = detect_checkpoint_variant(root / "upload/ParT_checkpoint.pt")
    assert variant == "legacy_mha"
    stored = np.fromfile(package / "data/model_weights.fp32.bin", dtype=np.float32)
    expected = np.concatenate(
        [
            np.ascontiguousarray(value.reshape(value.shape[0], -1), dtype=np.float32).reshape(-1)
            for name, value in state.items()
            if name.endswith(".weight") and value.ndim >= 2
        ]
    )
    assert np.array_equal(stored, expected)
    assert stored.tobytes() == expected.tobytes()

    expected_params = []
    consumed = set()
    for name, value in state.items():
        if not (name.endswith(".weight") and value.ndim >= 2):
            continue
        bias_name = name[:-7] + ".bias"
        bias = state.get(bias_name, np.zeros(value.shape[0], dtype=np.float32))
        expected_params.extend(np.ascontiguousarray(bias, dtype=np.float32).reshape(-1).tolist())
        consumed.add(name)
        if bias_name in state:
            consumed.add(bias_name)
    for name, value in state.items():
        if name in consumed or name.endswith("num_batches_tracked"):
            continue
        expected_params.extend(np.ascontiguousarray(value, dtype=np.float32).reshape(-1).tolist())
    stored_params = np.fromfile(package / "data/model_params.fp32.bin", dtype=np.float32)
    expected_params = np.ascontiguousarray(expected_params, dtype=np.float32)
    assert stored_params.tobytes() == expected_params.tobytes()

    rng = np.random.default_rng(20260827)
    features = rng.normal(size=(2, 15, 16)).astype(np.float32)
    vectors = rng.normal(size=(2, 4, 16)).astype(np.float32)
    vectors[:, 3] = np.sqrt((vectors[:, :3] ** 2).sum(1) + rng.uniform(.01, 2, (2, 16)))
    mask = np.zeros((2, 16), dtype=bool)
    mask[:, :8] = True
    logits = WeaverParT(state, model_variant=variant)(features, vectors, mask)
    assert np.isfinite(logits).all()

    printed_features = np.zeros((1, 15, 16), dtype=np.float32)
    printed_vectors = np.zeros((1, 4, 16), dtype=np.float32)
    printed_mask = np.zeros((1, 16), dtype=bool)
    printed_features[0, :, :14] = PRINTED_JET[:, 2:17].T
    printed_vectors[0, :, :14] = PRINTED_JET[:, 17:21].T
    printed_mask[0, :14] = True
    printed_probabilities = softmax(
        WeaverParT(state, model_variant=variant)(printed_features, printed_vectors, printed_mask)
    )
    np.testing.assert_allclose(
        printed_probabilities[0], [0.25344047, 0.74655953], rtol=0, atol=2e-6
    )

    old_checkpoint = root / "upload/net_best_epoch_state.pt"
    old_probabilities = None
    if old_checkpoint.exists():
        old_variant = detect_checkpoint_variant(old_checkpoint)
        assert old_variant == "modern_custom"
        old_probabilities = softmax(
            WeaverParT(
                load_state_dict(old_checkpoint), model_variant=old_variant
            )(printed_features, printed_vectors, printed_mask)
        )
        np.testing.assert_allclose(
            old_probabilities[0], [0.22407684, 0.77592320], rtol=0, atol=2e-6
        )
    print(
        {
            "weights": int(stored.size),
            "byte_exact": True,
            "weights_sha256": hashlib.sha256(stored.tobytes()).hexdigest(),
            "params_sha256": hashlib.sha256(stored_params.tobytes()).hexdigest(),
            "reference_probabilities": softmax(logits).tolist(),
            "printed_jet_probabilities": printed_probabilities.tolist(),
            "old_checkpoint_printed_jet_probabilities": (
                None if old_probabilities is None else old_probabilities.tolist()
            ),
        }
    )


if __name__ == "__main__":
    main()
