DY/W+jets prompt efficiency & fake rate for 3 IDs incl. modified HEEP Run3

Run:
  source env.sh
  python3 -u dyw_run3.py --max-files 100     # 0 = all files
  python3 dyw_run3.py --plots-only
Outputs: out/dyw_run3_<tag>/ (per-component .root, rates.csv, *_files.txt),
         imgs/dyw_run3_<tag>/*.{png,pdf}.

Self-contained Run3-dedicated script; does NOT modify the existing Study-2 scripts
(compare_prompt_heep.py, plot_dyw_cutflow.py, plot_id_inputs.py). Input production is
../26.10.08_22preEE (schema v6 ntuples that store ModifiedHEEPRun3ID); /pnfs 261008
directories are resolved from that folder's inputdataMC.dat. DY_Jets(+ext1/2) summed into
"DY"; W_Jets is "W". Unit electron weights; no xsec/lumi/pileup/signed MC weights.

Two SEPARATE measurements (not a comparison), EB and EE split:
  prompt efficiency = pass / baseline for reco electrons matched to a hard-process final-state
    gen electron (|eta_gen|<2.5) within dR<0.001.
  fake rate         = pass / baseline for reco electrons with NO eligible gen electron within
    dR>=0.1 (fake-candidate proxy, NOT a per-jet mis-ID; no W/Z ancestry stored).
Baseline: reco raw pT>20 GeV, EB |etaSC|<1.4442, EE 1.566<|etaSC|<2.5, no ID preselection
(existing minLeptonsSkim=2 remains). "pass" = ID 12-bit VID bitmap == 4095.

Three IDs per plot (same electrons): HEEP V70 (electron_heepV70Bitmap),
modified HEEP (electron_modHeepBitmap), modified HEEP Run3 (electron_bitmapModifiedHEEPRun3ID).
Per (group DY|W) x (population prompt|fake) x (region EB|EE):
  *_rate_vs_pt : passing fraction vs raw pT (Clopper-Pearson 68%), 3 ID curves.
  *_cutflow    : cumulative pass fraction through the 12 VID cuts, 3 ID curves.
rates.csv: integrated denominator/numerator/rate (+68% CP) at pT>20 and pT>35, per ID.

Note: in EB the modified and Run3 IDs are nearly identical (clean prompt barrel electrons);
the Run3 differences show up in EE, where it raises both prompt efficiency and fake rate.
Caveats as in prompt_heep_README.txt (not PF-isolated; HEEP V70 track iso not stored;
modEcalHcalIso ECAL-only).
