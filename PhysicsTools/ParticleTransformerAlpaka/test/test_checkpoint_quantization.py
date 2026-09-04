import pathlib,sys,numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).parents[1]/"scripts"))
from weaver_checkpoint import detect_checkpoint_variant,load_state_dict
from weaver_reference import WeaverParT,softmax

def main():
    root=pathlib.Path(__file__).parents[2]
    checkpoints=[root/"upload/net_best_epoch_state.pt",root/"upload/ParT_checkpoint.pt"]
    rows=np.loadtxt(pathlib.Path(__file__).with_name("printed_jet.txt"),usecols=range(1,22),dtype=np.float32)
    printed_features=np.zeros((1,15,16),np.float32);printed_vectors=np.zeros((1,4,16),np.float32);printed_mask=np.zeros((1,16),bool)
    printed_features[0,:,:14]=rows[:,2:17].T;printed_vectors[0,:,:14]=rows[:,17:21].T;printed_mask[0,:14]=True
    for checkpoint in checkpoints:
        if not checkpoint.exists():
            continue
        s=load_state_dict(checkpoint)
        variant=detect_checkpoint_variant(checkpoint)
        assert len(s)==235 and s["fc.0.weight"].shape==(2,128)
        assert "blocks.0.attn.in_proj.weight" in s
        assert "blocks.0.attn.in_proj_weight" not in s
        rng=np.random.default_rng(20260825); max_abs=[]; argmax_equal=[]
        fp=WeaverParT(s,False,variant); qi=WeaverParT(s,True,variant)
        for n in (1,4,8,16):
            x=rng.normal(size=(2,15,16)).astype(np.float32);v=rng.normal(size=(2,4,16)).astype(np.float32);v[:,3]=np.sqrt((v[:,:3]**2).sum(1)+rng.uniform(.01,2,(2,16)))
            mask=np.zeros((2,16),bool);mask[:,:n]=True;x[:,:,n:]=0;v[:,:,n:]=0
            a=fp(x,v,mask);b=qi(x,v,mask);max_abs.append(float(np.max(np.abs(a-b))));argmax_equal.extend(np.argmax(a,1)==np.argmax(b,1))
        printed_fp=softmax(fp(printed_features,printed_vectors,printed_mask));printed_q=softmax(qi(printed_features,printed_vectors,printed_mask))
        expected=0.25294158 if variant=="legacy_mha" else 0.22336726
        np.testing.assert_allclose(printed_q[0,0],expected,rtol=0,atol=2e-6)
        result={"checkpoint":checkpoint.name,"model_variant":variant,"max_abs_logit_error":max(max_abs),"mean_case_max_abs_error":float(np.mean(max_abs)),"argmax_agreement":float(np.mean(argmax_equal)),"printed_jet_fp32":printed_fp.tolist(),"printed_jet_int8":printed_q.tolist()}
        print(result)
        assert np.isfinite(max_abs).all()
if __name__=="__main__":main()
