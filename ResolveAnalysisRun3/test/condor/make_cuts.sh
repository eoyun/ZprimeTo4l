#!/bin/bash
# Per-cut marginal PASS efficiency for all 12 VID cuts, vs dR_gen and reco pT,
# EB/EE separate, 3 IDs (V70/modHEEP/modHEEP Run3), across all signal samples.
# -> imgs/cuts/<name>_<cutkey>_<region>[_pt].png
set -e
cd "$(dirname "$0")"
source env.sh
make plotCut
args=""
while read -r s; do
  [ -f "out/hists_${s}.root" ] && args="$args ${s}=out/hists_${s}.root" || echo "missing out/hists_${s}.root"
done < samples_signal.txt
[ -n "$args" ] || { echo "no histogram files found"; exit 1; }
CUTS="et sceta detainseed dphiin sigmaieta e2x5 hoe trkiso caloiso dxy misshits ecaldriven"
for cut in $CUTS; do
  for var in dr pt; do
    ./plotCut imgs/cuts "$cut" "$var" $args
  done
done