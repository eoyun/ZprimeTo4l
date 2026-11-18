#!/bin/bash
# Source inside the EL8 container; CMSSW supplies ROOT's dependency libraries.
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH="$MERGE_SCRAM_ARCH"
merge_previous_dir="$PWD"
cd "$MERGE_CMSSW_RELEASE/src" || return 1
merge_runtime=$(scram runtime -sh) || return 1
eval "$merge_runtime"
cd "$merge_previous_dir" || return 1
command -v hadd >/dev/null || return 1
