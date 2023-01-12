#!/bin/sh
export X509_USER_PROXY="{PROXYPATH}"
cd {CMSSWSETUP}
eval `scramv1 runtime -sh`
cd PhysicsTools/NanoAOD/test/{JOBFOLDER}
cmsRun {CFGNAME}
echo 'sending the file back'

python3 {DFNAME}

mkdir {OUTPUTPATH}
cp skimmed_{OUTPUTFILENAME} {OUTPUTPATH}
rm {OUTPUTFILENAME}

