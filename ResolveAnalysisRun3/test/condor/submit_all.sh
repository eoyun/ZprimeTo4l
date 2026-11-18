#!/bin/bash
# Build everything and submit one condor job per signal sample.
# Usage: ./submit_all.sh [PROD_DIR]   (default ../26.09.30_22preEE)
set -e
cd "$(dirname "$0")"
source env.sh
make
./build_filelists.sh "${1:-../26.10.08_22preEE}"
mkdir -p logs out
condor_submit submit.sub
echo "submitted. monitor: condor_q ; then plot with ./make_plots.sh"
