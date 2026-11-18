#!/bin/bash
# Build per-sample file lists of this production's CRAB output on /pnfs.
# Reads ../inputdata.dat + ../inputdataMC.dat, takes the task timestamp from each
# crab/<name>/crab_projects/crab_muonJpsi_<name>/.requestcache, and writes
# filelists/<name>.txt + samples.txt.
set -u
PROD=".."
BASE="/pnfs/knu.ac.kr/data/cms/store/user/yeo"
mkdir -p filelists logs ../hadd
: > samples.txt

cat "$PROD/inputdata.dat" "$PROD/inputdataMC.dat" | while read -r name dataset; do
  [ -z "$name" ] && continue
  primary=$(echo "$dataset" | cut -d/ -f2)
  cache="$PROD/crab/$name/crab_projects/crab_muonJpsi_$name/.requestcache"
  ts=$(grep -ao "[0-9]\{6\}_[0-9]\{6\}:yeo_crab_muonJpsi_$name" "$cache" 2>/dev/null | head -1 | cut -d: -f1)
  if [ -z "$ts" ]; then echo "WARN: no task timestamp for $name"; continue; fi
  dir="$BASE/$primary/muonJpsi_$name/$ts"
  find "$dir" -name '*.root' -not -path '*/log/*' 2>/dev/null | sort > "filelists/$name.txt"
  n=$(wc -l < "filelists/$name.txt")
  if [ "$n" -eq 0 ]; then echo "WARN: 0 files for $name ($dir)"; continue; fi
  printf "%-20s %5s files  (%s)\n" "$name" "$n" "$ts"
  echo "$name" >> samples.txt
done

echo "---"
echo "samples written: $(wc -l < samples.txt)  -> samples.txt"
