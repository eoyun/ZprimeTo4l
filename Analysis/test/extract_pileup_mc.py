#!/usr/bin/env python3
"""Extract unweighted BX=0 true-PU from unskimmed MC MiniAOD using CMSSW FWLite.

Run after cmsenv. Outputs a raw-count TH1D named pileup plus a normalized PNG.
This examines the supplied events; it does not establish their production profile.
"""

import argparse
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", nargs="+", help="MiniAOD ROOT paths or root:// URLs")
    inputs.add_argument("--file-list", type=Path, help="one MiniAOD path per line; # comments allowed")
    parser.add_argument("--output", type=Path, default=Path("data/pileup/pileup_mc_DY_2022preEE.root"))
    parser.add_argument("--max-events", type=int, default=-1, help="-1: all events; limit is across all input files")
    parser.add_argument("--bins", type=int, default=100)
    parser.add_argument("--max-pu", type=float, default=100.)
    parser.add_argument("--label", default="slimmedAddPileupInfo")
    parser.add_argument("--redirector", default="root://cms-xrd-global.cern.ch/", help="prefix for /store/ paths")
    args = parser.parse_args()
    if args.max_events == 0 or args.max_events < -1:
        parser.error("--max-events must be -1 or positive")
    if args.bins <= 0 or not math.isfinite(args.max_pu) or args.max_pu <= 0:
        parser.error("--bins and --max-pu must be positive and finite")
    files = args.input
    if args.file_list:
        files = [line.split("#", 1)[0].strip() for line in args.file_list.read_text().splitlines()]
        files = [line for line in files if line]
    if not files:
        parser.error("Empty input file list")
    files = [args.redirector.rstrip("/") + "/" + name if name.startswith("/store/") else name for name in files]
    output = args.output
    if output.suffix != ".root":
        parser.error("--output must end in .root")
    plot = output.with_suffix(".png")
    for path in (output, plot):
        if path.exists():
            parser.error("Output exists; choose another --output: {}".format(path))

    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    if ROOT.gSystem.Load("libFWCoreFWLite") < 0:
        raise RuntimeError("Cannot load FWLite; run inside CMSSW after cmsenv")
    ROOT.FWLiteEnabler.enable()
    ROOT.gSystem.Load("libDataFormatsFWLite")
    from DataFormats.FWLite import Events, Handle

    histogram = ROOT.TH1D("pileup", "MC true pileup;True mean interactions (BX=0);Events", args.bins, 0., args.max_pu)
    histogram.SetDirectory(0)
    histogram.Sumw2()
    handle = Handle("std::vector<PileupSummaryInfo>")
    events = Events(files, maxEvents=args.max_events)
    count = 0
    for event in events:
        if event.eventAuxiliary().isRealData():
            raise RuntimeError("Input contains real data; supply MC MiniAODSIM")
        event.getByLabel(args.label, handle)
        if not handle.isValid():
            raise RuntimeError("Missing {} at event {}".format(args.label, count))
        values = [float(pu.getTrueNumInteractions()) for pu in handle.product() if pu.getBunchCrossing() == 0]
        if len(values) != 1 or not math.isfinite(values[0]) or values[0] < 0:
            raise RuntimeError("Expected one valid BX=0 true-PU value at event {}: {}".format(count, values))
        # No selection, generator weight, PU weight, cross section or SF.
        histogram.Fill(values[0], 1.)
        count += 1
        if count % 10000 == 0:
            print("Processed {} events".format(count), flush=True)
    if count == 0:
        raise RuntimeError("No events read; no output produced")
    integral = histogram.Integral()
    if integral <= 0:
        raise RuntimeError("All events are outside the histogram range; increase --max-pu")
    summary = {
        "input_files": files, "label": args.label, "max_events": args.max_events,
        "events_read": count, "events_in_range": integral,
        "underflow": histogram.GetBinContent(0),
        "overflow": histogram.GetBinContent(args.bins + 1),
        "mean_in_range": histogram.GetMean(),
        "empty_bins_in_range": sum(histogram.GetBinContent(i) == 0 for i in range(1, args.bins + 1)),
        "weighting": "one count per event; no analysis cuts or weights",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    destination = ROOT.TFile.Open(str(output), "NEW")
    if not destination or destination.IsZombie():
        raise RuntimeError("Cannot create {}".format(output))
    destination.cd()
    histogram.Write()
    ROOT.TNamed("provenance", json.dumps(summary, sort_keys=True)).Write()
    destination.Close()
    normalized = histogram.Clone("pileup_normalized")
    normalized.SetDirectory(0)
    normalized.Scale(1. / integral)
    normalized.GetYaxis().SetTitle("Probability per bin (in-range normalized)")
    normalized.SetStats(False)
    canvas = ROOT.TCanvas("pileup_canvas", "MC pileup", 900, 700)
    normalized.Draw("E1 HIST")
    canvas.SaveAs(str(plot))
    print(json.dumps(summary, indent=2))
    if summary["overflow"] or summary["underflow"]:
        print("Nonzero flow bins: expand the range for both data and MC before making weights.")
    print("Created {} and {}".format(output.resolve(), plot.resolve()))


if __name__ == "__main__":
    main()
