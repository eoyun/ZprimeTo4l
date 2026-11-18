"""Check the ROOT schema produced by cmsRun ... maxEvents=0 (not event validation)."""
import sys
import ROOT

f = ROOT.TFile.Open(sys.argv[1])
assert f and not f.IsZombie(), "Cannot open output"
d = f.GetDirectory("mergedLeptonNtuplizer")
assert d, "Missing module directory"
tree = d.Get("Events")
assert tree and tree.GetEntries() == 0
branches = {b.GetName(): b.GetClassName() for b in tree.GetListOfBranches()}
for name in ("muon_pt", "muon_modifiedtrackiso", "muon_modifiedIsoRun2",
             "muon_emS25", "muon_segCompatibility", "electron_pt"):
    assert branches.get(name) == "vector<float>", (name, branches.get(name))
for name in ("muon_index", "muon_passMergedBaseline", "muon_passResolved",
             "electron_passRun2Pt", "electron_modifiedHeepRun3Valid", "electron_passModifiedHeepRun3"):
    assert branches.get(name) == "vector<int>", (name, branches.get(name))
for name in ("met_pt", "met_phi", "puppiMet_pt", "hlt_names", "filter_status",
             "trigObj_pathOffsets", "trigObj_labels", "gen_pt"):
    assert name in branches, name
assert tree.GetLeaf("event").GetTypeName() == "ULong64_t"
assert tree.GetLeaf("genWeight").GetTypeName() == "Double_t"
assert not any("bdt" in name.lower() or "m4l" in name for name in branches)
norm = d.Get("norm")
assert [norm.GetBinContent(i) for i in (1, 2, 3)] == [0., 0., 0.]
meta = d.Get("meta")
assert meta.GetBinContent(1) == 1
assert [meta.GetBinContent(i) for i in (3, 4, 5)] == [0., 0., 0.]
assert d.Get("configuration")
print("Empty ntuple schema check passed:", len(branches), "branches")
f.Close()
