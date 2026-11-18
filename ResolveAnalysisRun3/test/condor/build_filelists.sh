#!/bin/bash
# Build per-signal-sample file lists from the /pnfs ntuple output of a production folder.
# Usage: ./build_filelists.sh [PROD_DIR]   (default ../26.09.30_22preEE)
# Reads the HToAATo4L lines of PROD_DIR/inputdataMC.dat, resolves the /pnfs muonJpsi_<name>
# directory (latest timestamp), and writes filelists/<name>.txt + samples_signal.txt.
set -u
PROD="${1:-../26.10.08_22preEE}"
BASE="/pnfs/knu.ac.kr/data/cms/store/user/yeo"
DAT="$PROD/inputdataMC.dat"

[ -f "$DAT" ] || { echo "ERROR: $DAT not found"; exit 1; }
# Timestamp prefix from the production folder name, e.g. 26.10.08_22preEE -> 261008
TS=$(basename "$PROD" | grep -oE '[0-9]{2}\.[0-9]{2}\.[0-9]{2}' | tr -d '.')
[ -n "$TS" ] || { echo "ERROR: cannot derive timestamp from $PROD"; exit 1; }
echo "production timestamp prefix: ${TS}"
mkdir -p filelists
: > samples_signal.txt

grep -i HToAA "$DAT" | while read -r name dataset; do
  primary=$(echo "$dataset" | cut -d/ -f2)
  ts=$(ls -d "$BASE/$primary/muonJpsi_${name}"/${TS}* 2>/dev/null | sort | tail -1)
  if [ -z "$ts" ]; then echo "WARN: no /pnfs dir for $name ($primary)"; continue; fi
  find "$ts" -name 'Events*.root' 2>/dev/null | sort > "filelists/${name}.txt"
  n=$(wc -l < "filelists/${name}.txt")
  if [ "$n" -eq 0 ]; then echo "WARN: 0 files for $name"; continue; fi
  printf "%-32s %4s files  (%s)\n" "$name" "$n" "$(basename "$ts")"
  echo "$name" >> samples_signal.txt
done

echo "---"
echo "samples written: $(wc -l < samples_signal.txt)  -> samples_signal.txt"
