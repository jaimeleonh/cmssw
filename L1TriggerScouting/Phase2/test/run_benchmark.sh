minConstituents=2

#for model in deepset_cpp
#for model in softjet_part_unique
# for model in part_int8 deepset_int8 model_deepsets #softjet_deepsets
# for model in part_10k_12const
for model in part_2p5k part_10k deepsets_10k deepsets_2p5k  deepsets_50k   part_2p2m      part_300k      pnet_10k       pnet_2p5k      pnet_50k deepsets_2p2m  deepsets_300k  part_10k       part_2p5k      part_50k       pnet_2p2m      pnet_300k
do

    for backend in cuda_async # serial_sync # 
    do
        for numThreads in 1 #2 4 8 #16
        do
            for numFwkStreams in 1 #4 8 # 16
            do
                for batchSize in 128 256 512 1024 2048 4096 8192 16384 #32768 65536 #262144 524288 1048576 2097152 #1000000 2000000 4000000
                # for batchSize in 32 128 256 512 1024 2048 4096 8192 16384 65536 262144
                do
                    if [[ ! -f results_${model}/resources.cms-p2scout-ngt-01.${backend}.${numThreads}.${numFwkStreams}.${batchSize}.json ]]
                    then
    #                     tmp="ncu --kernel-name JetInferenceKernel \
    #     --metrics sm__pipe_fma_cycles_active.avg.pct_of_peak_sustained_active,\
    # sm__warps_active.avg.pct_of_peak_sustained_active,\
    # smsp__average_warps_issue_stalled_barrier_per_issue_active.ratio \
                        mkdir -p results_${model}
                        tmp="CUDA_VISIBLE_DEVICES=0 cmsRun runScoutingPhase2PF_cfg.py  buBaseDir=/mnt/ramdisk/$USER/raw fuBaseDir=/mnt/ramdisk/$USER/raw  buNumStreams=2 runNumber=42 maxEvents=100  run=sc4Alpaka backend=$backend   numThreads=$numThreads numFwkStreams=$numFwkStreams batchSize=$batchSize jetTagModel=${model} jetAlgo=SCNMSWeighted RSeed=0.4 RClu=0.4 minConstituents=${minConstituents} | tee log2026"
                        #tmp='LD_PRELOAD="$CMSSW_BASE/src/L1TriggerScouting/JetTagging/data/softjet_part_full_int8_ops/libscouting_int8.so${LD_PRELOAD:+:$LD_PRELOAD}" cmsRun runScoutingPhase2PF_cfg.py  buBaseDir=/mnt/ramdisk/$USER/raw fuBaseDir=/mnt/ramdisk/$USER/raw  buNumStreams=2 runNumber=41 maxEvents=100  run=sc4Alpaka backend=$backend   numThreads=$numThreads numFwkStreams=$numFwkStreams batchSize=$batchSize | tee log2026'
                        #tmp="cmsRun runScoutingPhase2PF_cfg.py  buBaseDir=/mnt/ramdisk/$USER/raw fuBaseDir=/mnt/ramdisk/$USER/raw  buNumStreams=1 runNumber=39 maxEvents=100  run=sc4Alpaka backend=$backend   numThreads=$numThreads numFwkStreams=$numFwkStreams batchSize=$batchSize | tee log2026"
                        echo $tmp
                        eval $tmp
                    fi
                done
            done
        done
    done
done


