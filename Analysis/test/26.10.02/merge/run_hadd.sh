#!/bin/bash
set -eo pipefail
if [[ $# -ne 1 ]]; then
    echo "Usage: run_hadd.sh job.json" >&2
    exit 2
fi
merge_script_dir=$(cd "$(dirname "$0")" && pwd)
if [[ ${MERGE_IN_CONTAINER:-0} != 1 ]]; then
    export MERGE_IN_CONTAINER=1
    exec /cvmfs/cms.cern.ch/common/cmssw-el8 -B /pnfs,/u,/tmp --command-to-run \
        bash "$merge_script_dir/run_hadd.sh" "$1"
fi
merge_settings=$(python3 - "$1" <<'PY'
import json, shlex, sys
j = json.load(open(sys.argv[1]))
for key, name in [('scram_arch', 'MERGE_SCRAM_ARCH'), ('cmssw_release', 'MERGE_CMSSW_RELEASE')]:
    print('export ' + name + '=' + shlex.quote(j[key]))
PY
)
eval "$merge_settings"
source "$merge_script_dir/setup_root.sh"
exec python3 "$merge_script_dir/check_root.py" "$1"
