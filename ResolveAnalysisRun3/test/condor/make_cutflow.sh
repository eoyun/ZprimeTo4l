#!/bin/bash
# Modified-HEEP cutflow (EB vs EE) across all signal samples, after condor jobs finish.
set -e
cd "$(dirname "$0")"
source env.sh
make plotCutflow
args=""
while read -r s; do
  [ -f "out/hists_${s}.root" ] && args="$args ${s}=out/hists_${s}.root" || echo "missing out/hists_${s}.root"
done < samples_signal.txt
[ -n "$args" ] || { echo "no histogram files found"; exit 1; }
mode="${1:-mod}"   # mod (modified HEEP) or v70 (standard HEEP V70)
./plotCutflow imgs/cutflow "$mode" $args
