#!/usr/bin/env python3
"""Convert data/MC true-pileup histograms to a correctionlib schema-v2 weight file.

Run inside CMSSW (PyROOT required). Supply data variations explicitly; missing
variations are not synthesized. Input histograms must have identical bin edges.
"""

import argparse
import json
import math
from pathlib import Path


def read_histogram(path, name):
    import ROOT

    source = ROOT.TFile.Open(str(path), "READ")
    if not source or source.IsZombie():
        raise ValueError("Cannot open ROOT file: {}".format(path))
    try:
        hist = source.Get(name)
        if not hist or not hist.InheritsFrom("TH1") or hist.GetDimension() != 1:
            raise ValueError("{}: missing one-dimensional histogram {}".format(path, name))
        size = hist.GetNbinsX()
        if hist.GetBinContent(0) != 0 or hist.GetBinContent(size + 1) != 0:
            raise ValueError("{}: nonzero underflow/overflow; fix the common bin range".format(path))
        edges = [float(hist.GetXaxis().GetBinLowEdge(i)) for i in range(1, size + 2)]
        values = [float(hist.GetBinContent(i)) for i in range(1, size + 1)]
        return edges, values
    finally:
        source.Close()


def normalize(values):
    if not values or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("Pileup distributions must contain finite, nonnegative values")
    total = sum(values)
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Pileup histogram has invalid integral")
    return [v / total for v in values]


def make_weights(data, mc):
    if len(data) != len(mc):
        raise ValueError("Data/MC bin counts differ")
    data, mc = normalize(data), normalize(mc)
    result = []
    for index, (d, m) in enumerate(zip(data, mc)):
        if m == 0 and d > 0:
            raise ValueError("MC is zero but data is positive in bin {}; inspect MC coverage".format(index + 1))
        value = d / m if m else 0.0
        if not math.isfinite(value):
            raise ValueError("Nonfinite PU weight in bin {}".format(index + 1))
        result.append(value)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--era", required=True, help="e.g. 2022preEE")
    parser.add_argument("--data-nominal", type=Path, required=True)
    parser.add_argument("--data-up", type=Path)
    parser.add_argument("--data-down", type=Path)
    parser.add_argument("--mc", type=Path, required=True, help="ROOT file with the uncorrected MC true-PU distribution")
    parser.add_argument("--data-hist", default="pileup")
    parser.add_argument("--mc-hist", default="pileup")
    parser.add_argument("--output", type=Path, help="default: data/pileup/puWeight<era>.json")
    args = parser.parse_args()
    if bool(args.data_up) != bool(args.data_down):
        parser.error("Provide both --data-up and --data-down, or neither")
    if not args.era or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in args.era):
        parser.error("--era must contain only letters, digits, underscores or hyphens")
    edges, mc = read_histogram(args.mc, args.mc_hist)
    content = []
    sources = {}
    for variation, path in [("nominal", args.data_nominal), ("up", args.data_up), ("down", args.data_down)]:
        if path is None:
            continue
        data_edges, data = read_histogram(path, args.data_hist)
        if len(edges) != len(data_edges) or any(not math.isclose(a, b, rel_tol=0, abs_tol=1e-9) for a, b in zip(edges, data_edges)):
            raise ValueError("{}: data/MC bin edges differ".format(path))
        sources[variation] = str(path.resolve())
        content.append({"key": variation, "value": {
            "nodetype": "binning", "input": "NumTrueInteractions", "flow": "clamp",
            "edges": edges, "content": make_weights(data, mc)}})
    variations = ", ".join(item["key"] for item in content)
    payload = {"schema_version": 2, "corrections": [{
        "name": "puWeight" + args.era, "version": 1,
        "description": "Normalized data / normalized MC true-PU ratio. Data: {}; MC: {} ({}). Both-zero bins map to zero; out-of-range inputs clamp to edge bins.".format(json.dumps(sources), args.mc.resolve(), args.mc_hist),
        "inputs": [
            {"name": "NumTrueInteractions", "type": "real", "description": "Number of true interactions (BX=0)"},
            {"name": "weights", "type": "string", "description": variations}],
        "output": {"name": "weight", "type": "real", "description": "Event weight for pileup reweighting"},
        "data": {"nodetype": "category", "input": "weights", "content": content}}]}
    encoded = json.dumps(payload, indent=2, allow_nan=False) + "\n"
    try:
        import correctionlib
    except ImportError:
        print("correctionlib unavailable: schema validation skipped")
    else:
        correctionlib.CorrectionSet.from_string(encoded)
    output = args.output or Path("data/pileup") / ("puWeight" + args.era + ".json")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        handle.write(encoded)
    print("Created {} (variations: {})".format(output.resolve(), variations))


if __name__ == "__main__":
    main()
