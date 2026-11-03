#!/usr/bin/env python3
"""Generate the synthetic jets and expected b-tag scores embedded in
runSoftJetTaggingONNXTest.py.

The jet constituents below are fixed by hand; this script derives from them:
  * the feature tensors, replicating step by step - in float32 - the arithmetic
    of ComputeParticleTransformerFeaturesKernelONNX
    (L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernelONNX.dev.cc),
  * the reference `softmax` scores, computed with the Python ONNX Runtime,
  * and prints them as a python block to paste into runSoftJetTaggingONNXTest.py.

Re-run (and update the cfg) whenever the model or the test data change:

  python3 generateSoftJetONNXTest.py /path/to/model.onnx
"""

import sys

import numpy as np
import onnxruntime as ort

F = np.float32


# ---------------------------------------------------------------------------
# synthetic event: 3 BX crossing, 5 jets
#   BX 0: jet0 (5 constituents), jet1 (20 -> truncated to 16 by the kernel), jet2 (4)
#   BX 1: jet3 (8 constituents), jet4 (1 constituent)
#   BX 2: empty (exercises the empty-BX code paths)
# one vertex per BX 0/1, none in BX 2
# ---------------------------------------------------------------------------

def make_constituents(pt, eta, phi, z0, dxy, puppiw, pdgid):
    return {
        "pt": np.array(pt, dtype=F),
        "eta": np.array(eta, dtype=F),
        "phi": np.array(phi, dtype=F),
        "z0": np.array(z0, dtype=F),
        "dxy": np.array(dxy, dtype=F),
        "puppiw": np.array(puppiw, dtype=F),
        "pdgid": np.array(pdgid, dtype=np.int32),
    }


jets = []

# jet0: 5 constituents
jets.append(make_constituents(
    pt=[45.0, 12.0, 3.5, 1.8, 0.9],
    eta=[0.11, 0.35, 1.9, -0.44, 2.2],
    phi=[0.05, 1.20, -0.60, 2.50, -1.70],
    z0=[0.10, 0.12, 0.35, -0.20, 0.05],
    dxy=[0.001, 0.010, 0.05, 0.20, 0.002],
    puppiw=[0.99, 0.85, 0.50, 0.30, 0.90],
    pdgid=[211, -211, 130, 22, 211],
))

# jet1: 20 constituents, the kernel keeps the first 16
rng = np.random.default_rng(42)
jets.append(make_constituents(
    pt=np.round(rng.uniform(0.6, 30.0, 20), 3).tolist(),
    eta=np.round(rng.uniform(-2.4, 2.4, 20), 3).tolist(),
    phi=np.round(rng.uniform(-3.0, 3.0, 20), 3).tolist(),
    z0=np.round(rng.uniform(-0.5, 0.5, 20), 3).tolist(),
    dxy=np.round(rng.uniform(0.0, 0.3, 20), 3).tolist(),
    puppiw=np.round(rng.uniform(0.1, 1.0, 20), 3).tolist(),
    pdgid=rng.choice([211, -211, 130, 22, 11, -13, 2112, 321, 2212], 20).tolist(),
))

# jet2: 4 constituents
jets.append(make_constituents(
    pt=[9.0, 4.0, 2.2, 1.1],
    eta=[-1.50, -1.20, -1.80, -0.90],
    phi=[2.80, 2.50, -2.90, 3.10],
    z0=[0.02, 0.08, -0.10, 0.15],
    dxy=[0.005, 0.020, 0.100, 0.050],
    puppiw=[0.95, 0.80, 0.60, 0.70],
    pdgid=[211, 22, 130, -211],
))

# jet3: 8 constituents (BX 1)
jets.append(make_constituents(
    pt=[25.0, 18.0, 6.5, 3.0, 2.0, 1.5, 1.0, 0.7],
    eta=[0.30, 0.10, 0.55, -0.20, 0.80, -0.60, 1.10, 0.00],
    phi=[-0.20, 0.15, -0.45, 0.70, -0.90, 1.30, -1.50, 0.05],
    z0=[-0.40, -0.38, -0.30, -0.50, -0.20, -0.45, 0.10, -0.42],
    dxy=[0.002, 0.004, 0.030, 0.150, 0.010, 0.060, 0.250, 0.001],
    puppiw=[1.00, 0.90, 0.70, 0.40, 0.85, 0.55, 0.20, 0.95],
    pdgid=[-13, 211, -211, 130, 22, 211, -211, 22],
))

# jet4: single constituent (BX 1)
jets.append(make_constituents(
    pt=[33.0],
    eta=[0.90],
    phi=[1.10],
    z0=[-0.35],
    dxy=[0.010],
    puppiw=[0.98],
    pdgid=[211],
))

jets_bx = [0, 0, 0, 1, 1]           # BX of each jet
vertex_z0 = [F(0.15), F(-0.40)]     # one vertex in BX 0 and BX 1
n_bx = 3
max_particles = 16


# ---------------------------------------------------------------------------
# replication of the kernel arithmetic (all in float32)
# ---------------------------------------------------------------------------

def clip(x, lo=-5.0, hi=5.0):
    return np.minimum(np.maximum(x, F(lo)), F(hi)).astype(F)


def safe_log(x):
    return np.log(np.maximum(x, F(1.0e-8))).astype(F)


def preprocess(v, sub, mul):
    return clip((v - F(sub)) * F(mul))


def jet_features(jet, vtx_z0):
    pt, eta, phi = jet["pt"], jet["eta"], jet["phi"]
    keep = pt > 0
    # keep only constituents with pt > 0, then at most max_particles of them
    idx = np.nonzero(keep)[0][:max_particles]
    pt, eta = pt[idx], eta[idx]
    phi, z0, dxy = jet["phi"][idx], jet["z0"][idx], jet["dxy"][idx]
    puppiw, pdgid = jet["puppiw"][idx], jet["pdgid"][idx]

    sin_phi = np.sin(phi).astype(F)
    cos_phi = np.cos(phi).astype(F)
    sinh_eta = np.sinh(eta).astype(F)
    cosh_eta = np.cosh(eta).astype(F)
    # pass 1: jet axis from the 4-momentum sum
    px = (pt * cos_phi).astype(F)
    py = (pt * sin_phi).astype(F)
    pz = (pt * sinh_eta).astype(F)
    e = (pt * cosh_eta).astype(F)
    ax_px, ax_py, ax_pz, ax_e = px.sum(), py.sum(), pz.sum(), e.sum()
    jet_pt = F(np.sqrt(F(ax_px * ax_px + ax_py * ax_py)))
    jet_eta = F(np.arcsinh(F(ax_pz / np.maximum(jet_pt, F(1.0e-8)))))
    jet_phi = np.arctan2(ax_py, ax_px).astype(F)

    n = len(idx)
    feats = np.zeros((15, n), dtype=F)
    vecs = np.zeros((4, n), dtype=F)
    mask = np.zeros((1, n), dtype=F)

    eta_sign = F(1.0) if jet_eta >= 0 else F(-1.0)
    deta = ((eta - jet_eta) * eta_sign).astype(F)
    dphi = (phi - jet_phi).astype(F)
    two_pi = F(2.0 * np.pi)
    dphi = np.where(dphi > F(np.pi), dphi - two_pi, dphi)
    dphi = np.where(dphi < F(-np.pi), dphi + two_pi, dphi)
    dr = np.sqrt((deta * deta + dphi * dphi).astype(F)).astype(F)

    feats[0] = preprocess(safe_log(pt), 1.7, 0.7)
    feats[1] = preprocess(safe_log(e), 2.0, 0.7)
    feats[2] = preprocess(safe_log((pt / np.maximum(jet_pt, F(1.0e-8))).astype(F)), -4.7, 0.7)
    feats[3] = preprocess(safe_log((e / np.maximum(ax_e, F(1.0e-8))).astype(F)), -4.7, 0.7)
    feats[4] = preprocess(dr, 0.2, 4.0)
    feats[5] = deta
    feats[6] = dphi
    feats[7] = (z0 - vtx_z0).astype(F)
    feats[8] = dxy
    feats[9] = puppiw
    abs_id = np.abs(pdgid)
    feats[10] = np.isin(abs_id, [211, 321, 2212]).astype(F)
    feats[11] = np.isin(abs_id, [130, 2112]).astype(F)
    feats[12] = (abs_id == 22).astype(F)
    feats[13] = (abs_id == 11).astype(F)
    feats[14] = (abs_id == 13).astype(F)

    vecs[0], vecs[1], vecs[2], vecs[3] = px, py, pz, e
    mask[0] = 1.0
    return feats, vecs, mask


def pad(a, n_front, n_back):
    return np.pad(a, [(0, 0), (n_front, n_front + n_back)])


def main(model_path):
    batch_feats, batch_vecs, batch_mask = [], [], []
    for jet, bx in zip(jets, jets_bx):
        feats, vecs, mask = jet_features(jet, vertex_z0[bx])
        batch_feats.append(pad(feats, 0, max_particles - feats.shape[1]))
        batch_vecs.append(pad(vecs, 0, max_particles - vecs.shape[1]))
        batch_mask.append(pad(mask, 0, max_particles - mask.shape[1]))

    # NOTE: this model cannot be run on a batch of jets: the pair_embed
    # Reshape node has a hardcoded [8, 136] target shape that only matches
    # batch size 1. Score therefore one jet at a time (as does the module
    # with maxBatchSize=1).
    session = ort.InferenceSession(model_path, providers=["CPUExecutionProvider"])
    scores = []
    for feats, vecs, mask in zip(batch_feats, batch_vecs, batch_mask):
        inputs = {
            "pf_features": feats[None].astype(np.float32),
            "pf_vectors": vecs[None].astype(np.float32),
            "pf_mask": mask[None].astype(np.float32),
        }
        scores.append(session.run(["softmax"], inputs)[0][0])
    scores = np.stack(scores)

    # flattened PF candidate list (association index is the identity)
    def flat(name):
        return [x for jet in jets for x in jet[name]]

    print("# --- begin generated block (generateSoftJetONNXTest.py) ---")
    for name in ["pt", "eta", "phi", "z0", "dxy", "puppiw"]:
        print(f"{name} = {flat(name)}")
    pdgids = flat("pdgid")
    print(f"pdgid = {[int(p) for p in pdgids]}")
    offsets = [0]
    for jet in jets:
        keep = np.nonzero(jet["pt"] > 0)[0].size
        offsets.append(offsets[-1] + keep)
    print(f"jetOffsets = {offsets}")
    bx_offsets = [0]
    for bx in range(n_bx):
        bx_offsets.append(bx_offsets[-1] + sum(b == bx for b in jets_bx))
    print(f"jetBxOffsets = {bx_offsets}")
    print(f"vertexZ0 = {[float(v) for v in vertex_z0]}")
    print(f"vertexBxOffsets = {[0, 1, 2, 2]}")
    print("expectedScores = [" + ", ".join(f"{s[0]:.9g}" for s in scores) + "]")
    print("# --- end generated block ---")
    print("\n# full softmax output (b, notB):", file=sys.stderr)
    for i, s in enumerate(scores):
        print(f"# jet{i}: {s[0]:.9g} {s[1]:.9g}", file=sys.stderr)


if __name__ == "__main__":
    model = sys.argv[1] if len(sys.argv) > 1 else \
        "L1TriggerScouting/JetTagging/data/btagging/softjet_part_example.onnx"
    import os
    for base in [os.environ.get("CMSSW_BASE", ""), os.environ.get("CMSSW_RELEASE_BASE", "")]:
        cand = os.path.join(base, "src", model) if base else None
        if cand and os.path.exists(cand):
            model = cand
            break
    main(model)
