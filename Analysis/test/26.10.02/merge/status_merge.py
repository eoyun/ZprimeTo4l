#!/usr/bin/env python3
"""Display sample-level merge progress from a prepared run's shared outputs."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run", type=Path)
args = parser.parse_args()
manifest = json.loads((args.run / "manifest.json").read_text())
done = ready = blocked = 0
for sample in manifest["samples"]:
    if sample["state"] != "ready":
        blocked += 1
        print(sample["sample"], "BLOCKED:", sample["reason"], sample.get("job_states", ""))
        continue
    ready += 1
    output = Path(sample["output"])
    if output.exists() and Path(str(output) + ".json").exists():
        done += 1
        print(sample["sample"], "VERIFIED:", output)
    else:
        print(sample["sample"], "PENDING (check DAG/Condor logs)")
print(f"Verified {done}/{ready}; blocked {blocked}")
