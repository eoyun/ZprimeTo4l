#!/bin/bash
# Local (no-condor) fallback: run the filler over every signal sample, then plot.
# Handy for quick turnaround since signal samples are small (<=50 files each).
# Usage: ./run_local.sh [PROD_DIR]   (default ../26.09.30_22preEE)
set -e
cd "$(dirname "$0")"
source env.sh
make
./build_filelists.sh "${1:-../26.10.08_22preEE}"
mkdir -p out
while read -r s; do
  ./fillHeepEffi "out/hists_${s}.root" "filelists/${s}.txt"
done < samples_signal.txt
./make_plots.sh
