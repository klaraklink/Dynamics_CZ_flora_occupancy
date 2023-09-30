#!/usr/bin/python3
import sys

with open("pladias_frequency_20230621.txt") as f:
    for line in f:
        name = line.strip()
        with open(f"jobs_mod/{name}.sh", "w", newline='\n') as of:
            contents = f"""#!/bin/bash
#PBS -N {name}
#PBS -l select=1:ncpus=2:mem=2gb:scratch_local=0b
#PBS -l walltime=96:00:00 
#PBS -j oe
#PBS -m n
# The 4 lines above are options for scheduling system: job will run 1 hour at maximum, 1 machine with 6 processors + 4gb RAM memory + 10gb scratch memory are requested
# Chybovy vystup pripoji ke standardnimu vystupu (-j oe) a posle mail pri zacatku, chybe, skonceni ulohy (-m abe).

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno2/home/klaraklink/pladias_mod_ht # substitute username and path to your real username and path

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR, specname = {name}" >> $DATADIR/jobs_info.txt

#loads application modules
module add r/4.0.0-gcc
export PATH=/storage/brno2/home/klaraklink/jags_install/bin:$PATH
export LD_LIBRARY_PATH=/storage/brno2/home/klaraklink/jags_install/lib:$LD_LIBRARY_PATH

cd $DATADIR

Rscript OCC_autojags_parKK2_plad_ht.R {name}

clean_scratch

"""
            of.write(contents + "\n")
            of.close()
