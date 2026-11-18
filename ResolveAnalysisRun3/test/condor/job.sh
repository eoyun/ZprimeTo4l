#!/bin/bash
# Condor executable. Runs the flat-ntuple HEEP-efficiency filler over one sample's files.
# transfer_input_files delivers: fillHeepEffi (binary) and <sample>.txt (file list).
set -e
sample="$1"
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc11-opt/setup.sh
chmod +x fillHeepEffi 2>/dev/null || true
./fillHeepEffi "hists_${sample}.root" "${sample}.txt"
