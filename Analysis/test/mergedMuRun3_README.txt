mergedMuRun3: Run-3 merged-muon analysis scaffold

Build in a CMSSW_13_0_13 runtime:
  scram b -j 4
Run (from this test directory, replace the input path):
  cmsRun runMergedMuRun3_cfg.py inputFiles=file:/path/to/miniAOD.root outputFile=mergedMuRun3.root maxEvents=100 isMC=True
Data: isMC=False (quality flags read from RECO instead of PAT).
HLT/filter switches: requireTrigger=False requireMETFilters=False for explicit diagnostics only.

Topology: two resolved muons + one merged-muon candidate. No additional
track, MET input, MET cut, missing momentum, or transverse-mass reconstruction.
METfilters are standard event-quality flags, not a MET selection.

Initial configurable selections:
- Good PV: non-fake, ndof>4, |z|<24 cm, transverse position<2 cm.
- Baseline: loose muon ID, muon pT>=20 GeV and |eta|<=2.4.
- Exactly three baseline muons by default (can relax to >=3).
- Resolved roles: HighPt OR TrackerHighPt ID, valid R03 track iso / muon pT<0.1.
  This initial selection uses raw isolation, not the old boosted subtraction.
- Merged role: baseline only. No high-pT ID or isolation cut.
- Event HLT OR; no offline trigger plateau or trigger-object matching yet.
- No electron veto, charge, resolved-pair dR/mass, or final signal-region cuts yet.
These are scaffold defaults, not finalized analysis selections.

Output:
- cutflow counts weighted EVENTS, never candidate rows.
- events has one row for each event passing trigger/quality/PV/multiplicity,
  including zero-candidate events. nCandidates counts all valid assignments.
- candidates has one row per distinct role assignment. Pair order is by pT.
  Original slimmedMuons indices, run/lumi/event link roles and event rows.
  Multiple rows per event are intentional: do not treat them as independent
  event yields or sum their full event weights as a final physics distribution.
- bdtScore=-999, bdtValid=false. No fake score or automatic ID pass is assigned.
- Track iso uses isolationR03().sumPt in GeV; invalid values are -999.
- Merged calorimeter energy fields are saved as initial inspection inputs.
  The final trained feature list and preprocessing remain to be implemented.
- p4 values are unmodified pat::Muon p4. visibleMass is the three-object visible
  invariant mass, not a calibrated four-muon resonance mass. A merged momentum
  prescription must be established separately from merged identification.
- MC weight is the full generator weight; data weight is 1. PU, muon SF,
  momentum corrections, FF, luminosity normalization and systematics are pending.
  No Run-2 payload is silently applied.

Next implementation points:
1. Replace mergedId() with inference and exact training feature preprocessing.
2. Establish merged-ID working point and event-level candidate arbitration.
3. Validate baseline/overlap/isolation, triggers, merged p4 and event selections.
4. Add Run-3 corrections, efficiencies, background strategy and systematic shifts.
