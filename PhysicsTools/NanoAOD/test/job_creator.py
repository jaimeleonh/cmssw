import sys, os
from subprocess import call
cfg_name = "NANO_NANO_template.py"
executable_name = "executable_template.sh"
df_name = "df_template.py"
user_proxy = "/afs/cern.ch/work/j/jleonhol/private/L1/nanoaod/CMSSW_12_6_0_patch1/src/x509up_u124088"
cmssw_setup = os.path.join(os.environ["CMSSW_BASE"], "src/")


import argparse
parser = argparse.ArgumentParser(description='options')
parser.add_argument('-i','--input', dest='input', required=True)
parser.add_argument('-o','--output', dest='output', required=True)
parser.add_argument('-n','--nfiles', type=int, dest='nfiles', default=1)
options = parser.parse_args()

input_files_filename = options.input
output_folder = options.output

with open(input_files_filename) as f:
    input_files = [filename.strip() for filename in f.readlines()]
job_folder = input_files_filename.split(".")[0]

with open(cfg_name) as f:
    cfg = f.read()

with open(df_name) as f:
    df = f.read()

with open(executable_name) as f:
    executable = f.read()

if os.path.isdir(job_folder):
    print("Remove job folder %s before running" % job_folder)
    sys.exit()

os.mkdir(job_folder)

# for ifile, f in enumerate(input_files):
for ifile in range(0, len(input_files), options.nfiles):
    if options.nfiles == 1:
        f = '"%s"' % input_files[ifile]
    else:
        max_ifile = min(ifile + options.nfiles, len(input_files))
        f = ", ".join(['"%s"' % elem for elem in input_files[ifile:max_ifile]])
        ifile = ifile // options.nfiles
    this_nano = cfg.replace("{INPUTFILENAME}", f).replace("{OUTPUTFILENAME}", "ntuple_%s.root" % ifile)
    os.mkdir("%s/Job_%s/" % (job_folder, ifile))
    with open("%s/Job_%s/cfg.py" % (job_folder, ifile), "w+") as f2:
        f2.write(this_nano)
    this_df = df.replace("{OUTPUTFILENAME}", "ntuple_%s.root" % ifile
        ).replace("{OUTPUTPATH}", output_folder)
    with open("%s/Job_%s/df.py" % (job_folder, ifile), "w+") as f2:
        f2.write(this_df)
    this_executable = executable.replace("{PROXYPATH}", user_proxy
        ).replace("{CMSSWSETUP}", cmssw_setup
        ).replace("{JOBFOLDER}", "%s/Job_%s/" % (job_folder, ifile)
        ).replace("{CFGNAME}", "cfg.py"
        ).replace("{OUTPUTFILENAME}", "ntuple_%s.root" % ifile
        ).replace("{OUTPUTPATH}", output_folder
        ).replace("{DFNAME}", "df.py")
    with open("%s/Job_%s/executable.sh" % (job_folder, ifile), "w+") as f2:
        f2.write(this_executable)
    
