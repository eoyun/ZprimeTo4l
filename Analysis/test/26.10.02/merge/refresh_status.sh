#!/bin/bash
# Prepare a DAG using fresh CRAB job states (requires an existing valid proxy).
set -eo pipefail
merge_here=$(cd "$(dirname "$0")" && pwd)
merge_config=${MERGE_CONFIG:-"$merge_here/config.json"}
merge_settings=$(python3 - "$merge_config" <<'PY'
import json, shlex, sys
c = json.load(open(sys.argv[1]))
for key, name in [('scram_arch', 'MERGE_SCRAM_ARCH'), ('cmssw_release', 'MERGE_CMSSW_RELEASE')]:
    print('export ' + name + '=' + shlex.quote(c[key]))
PY
)
eval "$merge_settings"
source "$merge_here/setup_root.sh"
source /cvmfs/cms.cern.ch/crab3/crab.sh prod
exec python3 "$merge_here/prepare_merge.py" --config "$merge_config" --refresh-status "$@"
