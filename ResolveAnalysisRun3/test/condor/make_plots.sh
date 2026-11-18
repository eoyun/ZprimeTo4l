#!/bin/bash
# Overlay modified-HEEP efficiency vs dR_gen across all signal samples (after condor jobs finish).
# Reads out/hists_<sample>.root for every sample in samples_signal.txt.
set -e
cd "$(dirname "$0")"
source env.sh
make plotHeepEffi
args=""
while read -r s; do
  [ -f "out/hists_${s}.root" ] && args="$args ${s}=out/hists_${s}.root" || echo "missing out/hists_${s}.root"
done < samples_signal.txt
[ -n "$args" ] || { echo "no histogram files found"; exit 1; }
mode="${1:-dr}"   # dr (efficiency vs dR_gen) or pt (efficiency vs reco pT)
./plotHeepEffi imgs/signal_all "$mode" $args
