#!/bin/bash
# Condor executable: hadd all files of one sample (list delivered as <sample>.txt).
set -e
sample="$1"
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc11-opt/setup.sh
echo "host: $(hostname)  files: $(wc -l < ${sample}.txt)"
hadd -f -j "${NCPU:-4}" "${sample}.root" $(cat "${sample}.txt")
ls -l "${sample}.root"
