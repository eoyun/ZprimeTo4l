# MergedAnalysisRun3 lepton ntuplizer

MiniAOD -> `mergedLeptonNtuplizer/Events`, one row per retained event and flat
`muon_*` / `electron_*` arrays. This package only produces inputs for downstream
merged ID and 2l+MM analysis. BDT inference, candidate building, missing-muon p4,
m3l/m4l, CR/SR selection, SF and momentum corrections are not implemented here.

## Selection contract (schema v1)

- PV: first `offlineSlimmedPrimaryVertices` entry, as in Run-2 `Analysis`.
  Events with no PV are dropped; no extra good-PV cut is introduced.
- Muon storage is the UNION of the two roles below. Original MiniAOD index is
  retained; global HighPt and TrackerHighPt flags are stored separately.
  - `passMergedBaseline`: global HighPt, reco pT >= 20 GeV, |reco eta| <= 2.4,
    following `Analysis/plugins/MergedMuCRanalyzer.cc`. No isolation cut.
  - `passResolvedAcceptance`: HighPt OR TrackerHighPt, TuneP pT >= 20 GeV,
    |reco eta| <= 2.4, following `Analysis/plugins/ResolvedMuCRanalyzer.cc`.
  - `passResolved`: resolved acceptance plus valid Run-2 modified track
    isolation / TuneP pT < 0.1. This is an object flag, not an event/pair selection.
  - Objects failing both HighPt IDs are not stored. This production does NOT
    support a later loose-muon baseline or Run-2 nonHighPt fake denominator.
- Electron storage: valid supercluster, |etaSC| <= 2.5, excluding the OPEN
  interval 1.4442 < |etaSC| < 1.566 (Run-2 pre-ID acceptance). All pT are retained;
  `passRun2Pt` records pT > 20 GeV (the modified-HEEP MinPtCut). No electron ID is applied. In particular,
  `modifiedHeepRun3Valid=0`, `passModifiedHeepRun3=-1` mean NOT EVALUATED.
  These are not selected-electron flags. Run-3 modified HEEP producers and
  their ValueMaps will be integrated after the study; current basic features
  alone are not sufficient to reproduce that ID offline.
- No HLT, event-filter, MET, mass, charge, isolation, channel or multiplicity
  hard cut (except optional `minLeptons`, default 0). Electron-only events remain.
  The ntuple therefore preserves the proposed low-mass CR without importing
  the old m3mu >= 50 GeV / MET > 50 GeV event cuts.

## Muon quantities

`muon_pt/eta/phi/mass/energy` are the unmodified PAT/reco p4; `muon_tuneP_*`
are separate. Do not treat one as the other. Inner/outer/global/best track
kinematics and validity, track-quality/ID inputs, calorimeter energy,
segment/calo compatibility and kink variables follow `MergedLepton/plugins/MergedMuon.cc`.
No trained feature list or model is required to produce this ntuple.

Two distinct modified isolations are stored in GeV:

| Branch | Companion pool |
| --- | --- |
| `muon_modifiedtrackiso` | All input muons with valid inner/TuneP tracks, as in MergedMuon/MuonJpsi studies |
| `muon_modifiedIsoRun2` | Resolved acceptance + HighPt/TrackerHighPt pool, as in Run-2 resolved selection |

Both subtract only the inner-track pT of the eligible companion with highest
TuneP pT. The matching uses reco muon vs companion inner-track directions,
inclusive 0.01 <= dR <= 0.3, |delta vz| <= 0.2, and SIGNED dxyBS <= 0.1.
Negative isolation is retained. Companion source indices are stored; a study
companion can be absent from the stored muon array. Compute before filtering.
Missing beam spot or invalid raw isolation gives -999 and `modifiedIsoValid=0`;
invalid isolation can never pass the resolved flag. Equal-pT ties use input order.

Optional track quantities use -999 (float) / -1 (int). Calorimeter and combined
quality values are kept as returned by MiniAOD, matching the study; their
validity flags are stored separately. No invalid BDT score is fabricated.

## Event information

- `genWeight` is a double, with full signed generator weight (data=1).
  `norm` stores Nevents, sumGenWeights, sumGenWeights2 BEFORE PV/skim rejection.
- `met_*` / `puppiMet_*` store pT, ORIGINAL phi and sumEt from `slimmedMETs` /
  `slimmedMETsPuppi`. They are the embedded nominal values; no correction or
  missing-muon substitution is applied. Missing/empty MET has valid=false and
  -999 values. MET systematic shifts are not stored in this first version.
- Full HLT and event-filter names/status arrays: -1=not run, 0=reject, 1=accept,
  2=error. Do not interpret nonzero as accepted. No menu/era choice is imposed.
- Trigger objects retain p4 direction, PDG id, trigger types, filter labels and paths satisfying
  `hasPathName(path,true,true)`. `typeOffsets`, `pathOffsets` and `labelOffsets` have Nobjects+1
  entries; object i uses flat entries [offset[i], offset[i+1]). Match only to
  paths accepted in `hlt_status`. Required trigger products fail loudly if absent.
- MC hard-process final-state gen electrons/muons are stored independently of
  reco selection for downstream matching. No truth label selects reco objects.
  Data does not read generator, pileup or gen particle products.
- `meta` records schema v1, isMC, and unapplied correction/ID flags.
  `configuration` is a TNamed containing the complete analyzer ParameterSet.

## Build and run

From a CMSSW_13_0_13 runtime in this working area:

```bash
scram b -j 4 ZprimeTo4l/MergedAnalysisRun3 testMergedRun3MuonIsolation
edmPluginRefresh "$CMSSW_BASE/lib/$SCRAM_ARCH/pluginZprimeTo4lMergedAnalysisRun3Plugins.so"
cmsRun ZprimeTo4l/MergedAnalysisRun3/test/runLeptonNtuplizer_cfg.py maxEvents=0 outputFile=/tmp/leptons_empty.root
python3 ZprimeTo4l/MergedAnalysisRun3/test/checkEmptyNtuple.py /tmp/leptons_empty.root
"$CMSSW_BASE/test/$SCRAM_ARCH/testMergedRun3MuonIsolation"
```

Paths above assume execution from `$CMSSW_BASE/src`. For real input:

```bash
cmsRun ZprimeTo4l/MergedAnalysisRun3/test/runLeptonNtuplizer_cfg.py inputFiles=file:/path/to/MiniAOD.root isMC=True maxEvents=100 outputFile=Leptons.root
```

For data use `isMC=False`. Event-filter process defaults to PAT for MC and RECO
for data; override `filterProcess=PAT` if the actual production requires it.
Only embedded quantities are read: no GlobalTag, conditions-dependent producer,
HEEP producer, BDT model, or Run-2 correction payload is needed.

The zero-event smoke check validates loading, configuration and ROOT booking;
it does not validate MiniAOD event reading or physics yields. A representative
data and MC event sample will still be needed before production.
