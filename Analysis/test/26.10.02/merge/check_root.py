#!/usr/bin/env python3
"""Run hadd and verify tree entries and normalization histograms before publication."""
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import ROOT

ROOT.gROOT.SetBatch(True)
NORMS = {"totWeightedSum", "sumWeights", "sumGenWeights", "nEvents", "cutflow"}


def inspect(path):
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"Invalid or recovered ROOT file: {path}")
    trees, norms, schema = {}, {}, {}

    def visit(directory, prefix=""):
        for key in directory.GetListOfKeys():
            name = prefix + key.GetName()
            if name in schema:
                continue  # ROOT key cycles refer to the same logical object.
            schema[name] = key.GetClassName()
            obj = directory.Get(key.GetName())
            if obj.InheritsFrom("TDirectory"):
                visit(obj, name + "/")
            elif obj.InheritsFrom("TTree"):
                trees[name] = int(obj.GetEntries())
                schema[name] = (key.GetClassName(), tuple(
                    (b.GetName(), b.GetTitle(), b.GetClassName())
                    for b in obj.GetListOfBranches()))
            elif obj.InheritsFrom("TH1") and key.GetName() in NORMS:
                norms[name] = [obj.GetBinContent(i) for i in range(obj.GetNcells())]
    try:
        visit(f)
        if not trees and not norms:
            raise RuntimeError(f"No analysis trees or normalization histograms: {path}")
        return schema, trees, norms
    finally:
        f.Close()


def main():
    job_path = Path(sys.argv[1]).resolve()
    job = json.loads(job_path.read_text())
    inputs = Path(job["input_list"]).read_text().splitlines()
    if not inputs or len(inputs) != len(set(inputs)):
        raise RuntimeError("Empty or duplicate input list")
    output = Path(job["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive reservation prevents duplicate workers from publishing concurrently.
    lock = Path(str(output) + ".lock")
    with lock.open("x"):
        pass
    try:
        if output.exists():
            raise RuntimeError(f"Output already exists; refusing overwrite: {output}")
        schema = None
        expected_trees, expected_norms = {}, {}
        for path in inputs:
            current, trees, norms = inspect(path)
            if schema is None:
                schema = current
                expected_trees = dict(trees)
                expected_norms = {k: list(v) for k, v in norms.items()}
            else:
                if current != schema:
                    raise RuntimeError(f"ROOT object or branch schema differs: {path}")
                for key, value in trees.items():
                    expected_trees[key] += value
                for key, values in norms.items():
                    if len(values) != len(expected_norms[key]):
                        raise RuntimeError(f"Normalization binning differs: {path}: {key}")
                    expected_norms[key] = [a + b for a, b in zip(expected_norms[key], values)]
        scratch = os.environ.get("_CONDOR_SCRATCH_DIR", os.environ.get("TMPDIR", "/tmp"))
        with tempfile.TemporaryDirectory(prefix="root-merge-", dir=scratch) as temp:
            target = Path(temp) / "merged.root"
            cmd = ["hadd", "-n", "50", str(target), "@" + job["input_list"]]
            print("Running:", cmd, flush=True)
            subprocess.run(cmd, check=True)
            actual_schema, trees, norms = inspect(target)
            if actual_schema != schema or trees != expected_trees:
                raise RuntimeError("Merged object schema or tree entry count differs")
            for key, values in expected_norms.items():
                if any(not math.isclose(a, b, rel_tol=1e-5, abs_tol=1e-5)
                       for a, b in zip(values, norms[key])):
                    raise RuntimeError(f"Merged normalization differs: {key}")
            # Stage on destination filesystem, then atomically publish without overwrite.
            with tempfile.NamedTemporaryFile(dir=output.parent, prefix=".merge-", delete=False) as staged:
                stage = Path(staged.name)
            try:
                import shutil
                shutil.copyfile(target, stage)
                os.link(stage, output)
            finally:
                stage.unlink(missing_ok=True)
            report = {"output": str(output), "inputs": inputs, "tree_entries": trees,
                      "normalization": norms, "root_version": ROOT.gROOT.GetVersion()}
            Path(str(output) + ".json").write_text(json.dumps(report, indent=2) + "\n")
            print("Verified and published:", output, flush=True)
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
