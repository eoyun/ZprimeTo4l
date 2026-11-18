# Modified-HEEP electron-ID studies (Run-3, flat ntuple, condor)

Standalone flat-ntuple analysis that runs natively on el9 condor (no CMSSW container):
builds compiled C++ against an LCG view and reads schema-v6 `Events.root` ntuples
(`resolvedNtuplizer/Events`, with `gen_*` branches and the `ModifiedHEEPRun3ID` VID result)
directly from the locally-mounted /pnfs. Current input production: `../26.10.08_22preEE`.

## Environment
```
source env.sh      # /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc11-opt  (ROOT 6.28, gcc11)
make               # builds fillHeepEffi, plotHeepEffi, plotCutflow, runHeepEffi
```
Ntuples: `/pnfs/knu.ac.kr/data/cms/store/user/yeo/<primaryDataset>/muonJpsi_<name>/<timestamp>/000X/Events_*.root`.

## Run
```
./submit_all.sh ../26.10.08_22preEE   # build filelists (signal) + submit 1 condor job/sample
./make_plots.sh dr                    # efficiency vs dR_gen  (EB/EE split; V70 + modHEEP + modHEEP Run3)
./make_plots.sh pt                    # efficiency vs reco pT
./make_cutflow.sh mod                 # modified-HEEP sequential cutflow (EB vs EE)
./make_cutflow.sh v70                 # standard HEEP V70 cutflow
./make_cutflow.sh modrun3             # modified-HEEP Run3 cutflow
./run_local.sh                        # no-condor fallback (small samples)
```
(`build_filelists.sh` derives the /pnfs timestamp prefix from the production folder name,
e.g. `26.10.08_22preEE` -> `261008`.)
Outputs: `out/hists_<sample>.root` (per-sample histograms), `imgs/signal_all/` (efficiency),
`imgs/cutflow/` (cutflows).

## Code map
| file | role |
|---|---|
| `HeepEffiCore.h` | single source of truth: selection, dR/pT binning, cutflow fill |
| `fillHeepEffi.cc` | condor worker: filelist -> histograms |
| `plotHeepEffi.cc` | efficiency plots, mode `dr`\|`pt`, EB/EE separate; 3 curves (V70, modHEEP, modHEEP Run3) |
| `plotCutflow.cc` | sequential cutflow, mode `mod`\|`v70`\|`modrun3`, EB vs EE |
| `build_filelists.sh` | resolve /pnfs dirs from a production folder's inputdataMC.dat |
| `submit.sub`,`job.sh`,`submit_all.sh` | condor (1 job/sample) |
| `make_plots.sh`,`make_cutflow.sh`,`run_local.sh` | orchestration |
| `runHeepEffi.cc` | single-file before/after demo |

---

## Study 1 — Signal modified-HEEP ID efficiency  (DONE)

Signal H->AA->4l, A->ee, 15 mass points (HToAATo4L_H{250,750,2000}A{0p4,1,2,5,10}_22).
Measures **ID efficiency conditional on reconstruction**, EB and EE separately.

Selection (`HeepEffiCore.h::processTree`), per the AN2023_143 §4.2 baseline:
- **Event level:** none (no HLT / MET-filter / PV / Z-mass / trigger-matching / tag&probe).
- **Gen pair:** gen e+/e- paired by greedy nearest opposite-charge dR; require both
  |eta_gen|<2.5 and leading pT_gen>20 GeV; pair carries dR_gen.
- **Denominator (reco e):** ptRaw>10, matched to a selected-pair gen e with dR(reco,gen)<0.001
  (closest, each reco counted once); region by |etaSC|: EB<1.4442, EE 1.566-2.5 (gap excluded);
  no ID applied.
- **Numerator (3 IDs):** `electron_passHeepV70` (nominal V70),
  `electron_passModHeep` (modified HEEP), `electron_passModifiedHEEPRun3ID` (modified HEEP Run3).
  Unweighted. Cutflows use the matching 12-bit maps (`*_heepV70Bitmap`, `*_modHeepBitmap`,
  `*_bitmapModifiedHEEPRun3ID`).

Observables (each: 15 samples + 3 H-mass summed + 1 all; EB/EE separate; log-x):
- efficiency vs **dR_gen** (variable binning 1e-4..1.0) — 3 curves per canvas
- efficiency vs **reco pT** (10..1000 GeV) — 3 curves per canvas
- **cutflow** (survival fraction through the 12 VID cuts in order), EB vs EE, one set each for
  modified (`mod`), V70 (`v70`), and Run3 (`modrun3`).

Key results:
- Modified HEEP >> nominal V70 (all-signal EB ~0.84 vs ~0.32).
- **EE:** nominal V70 collapses in mid-dR_gen (down to ~0.1); the old modified HEEP only
  partially recovers (~0.35-0.6); **modified HEEP Run3 recovers strongly (~0.77-0.86)**.
  So the Run3 version fixes the endcap deficit the old modified HEEP left. In EB the two
  modified IDs are effectively identical (clean prompt barrel electrons already recovered).
- The old EE loss traces (via cutflow) to dEtaInSeed (2nd-track recovery was EB-only) and the
  EE-only sigmaIEtaIEta<0.03; the Run3 producers address the endcap case.
- Modified gains also come from ET>20 (vs V70 35), E2x5/E5x5 bypass with a 2nd track, and
  modified track/ECAL isolation.
- Caveat: V70 uses ET>35 while both modified IDs use ET>20 — see AN §6.

---

## Study 2 — DY / W+jets: fake rate & prompt-electron ID  (DONE, python)

Standalone PyROOT/RDataFrame scripts (run locally, not the C++ condor path). Reads the stored
VID bitmaps (`electron_heepV70Bitmap` / `electron_modHeepBitmap`), cross-checked vs pass flags
and 4095. **Two separate measurements, NOT a nominal-vs-modified comparison.**

- `compare_prompt_heep.py` — engine. `--population matched|unmatched`, `--match-dr`, `--max-files`.
  DY_Jets(+ext1/2), W_Jets; 100 evenly-spaced files/production (or all); DY summed, unit weights.
  **prompt(matched):** reco e matched to a hard-process final-state gen e within dR<0.001.
  **fake(unmatched):** no eligible gen e within dR (>=0.1). Baseline reco pT>20 (and 35 for a
  common-ID comparison), EB/EE by etaSC, no ID preselection. Produces cut impact
  (individual / cumulative / N-1), input-variable distributions, matching same-subset, cut
  transitions, passing-rate vs pT (Clopper-Pearson), + CSVs (summary/cuts/flows/rates).
  -> out|imgs/{prompt_heep_100files, heep_unmatched_dr0p1_100files}.
  See `prompt_heep_README.txt`, `unmatched_heep_README.txt`.
- `plot_dyw_cutflow.py` — cumulative VID cutflow only; prompt+unmatched x EB/EE x pT{20,35},
  nominal (blue dashed) vs modified (red solid), baseline normalized to 1; reuses the same file
  manifests and validates counts against compare_prompt_heep. -> imgs/dyw_cutflow_only.
  See `dyw_cutflow_README.txt`.
- `plot_id_inputs.py` — SIGNAL ID-input distributions (dPerpIn, dEtaInSeed, sigmaIEtaIEta; signed;
  EB/EE) using the HeepEffiCore denominator; per-sample pages + compare_H*/compare_A* overlays.
  -> imgs/id_inputs. See `id_inputs_README.txt`.
- `dyw_run3.py` — Run3-dedicated (self-contained; does not touch the scripts above). DY/W prompt
  efficiency & fake rate for 3 IDs (V70, modified HEEP, modified HEEP Run3), EB/EE split, on the
  26.10.08 (schema-v6) production. Per group x population x region: rate-vs-pT (Clopper-Pearson)
  and cumulative cutflow, + rates.csv. -> out|imgs/dyw_run3_<tag>. See `dyw_run3_README.txt`.
  Finding: EB modified≈Run3; in EE Run3 raises both prompt efficiency and fake rate.

Caveats (documented in the READMEs): not PF-isolated; ntuple minLeptonsSkim=2 biases the
population (esp. W); no W/Z ancestry stored (prompt = hard-process-matched proxy, unmatched =
fake proxy, not a per-jet mis-ID); HEEP V70 track iso (dr03TkSumPtHEEP) is not stored, so
isolation effects come from the VID bits; modEcalHcalIso is ECAL-only.
