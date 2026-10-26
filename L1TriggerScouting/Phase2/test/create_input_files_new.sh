export APPTAINER_CACHEDIR=/tmp/$USER/cache
IMAGE=docker://gitlab-registry.cern.ch/scouting-demonstrator/phase-2/simple-clients/el9:3.2
## PF candidates from ZeroBias
# mkdir -p /mnt/ramdisk/$USER/raw/run000038 || true
# apptainer exec -B /eos -B /mnt/ramdisk $IMAGE p2scout-generator CMSSW /data/PF_Barrel.dump /mnt/ramdisk/$USER/raw/run000038/run000038_ls0001_index000001_stream00.raw --orbits 300 -T 1 -O 16384
# apptainer exec -B /eos -B /mnt/ramdisk $IMAGE p2scout-generator CMSSW /data/PF_Endcap.dump /mnt/ramdisk/$USER/raw/run000038/run000038_ls0001_index000001_stream01.raw --orbits 300 -T 1 -O 16384
# ## Puppi candidates from TTbar
# mkdir -p /mnt/ramdisk/$USER/raw/run000039 || true
# apptainer exec -B /eos -B /mnt/ramdisk $IMAGE  p2scout-generator CMSSW /eos/cms/store/cmst3/group/l1tr/gpetrucc/l1scout/binaries/puppi_TTbar_PU200.125X_v1.0.dump /mnt/ramdisk/$USER/raw/run000039/run000039_ls0001_index000001_stream00.raw --orbits 1000 -T 1
mkdir -p /mnt/ramdisk/$USER/raw/run000041 || true
apptainer exec -B /eos -B /mnt/ramdisk $IMAGE  p2scout-generator CMSSW puppi_mb.dump /mnt/ramdisk/$USER/raw/run000042/run000042_ls0001_index000001_stream00.raw --orbits 1000 -T 1 -s 0
apptainer exec -B /eos -B /mnt/ramdisk $IMAGE  p2scout-generator CMSSW vertices_mb.dump /mnt/ramdisk/$USER/raw/run000042/run000042_ls0001_index000001_stream01.raw --orbits 1000 -T 1 -s 0
