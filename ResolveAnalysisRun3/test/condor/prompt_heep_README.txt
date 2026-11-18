DY+jets / W+jets nominal and modified HEEP diagnostic

Run from this directory:
  source env.sh
  python3 -u compare_prompt_heep.py --max-files 100
  python3 compare_prompt_heep.py --max-files 100 --plots-only
Use --max-files 0 for all files (separate output directory).

Input filelists: DY_Jets, DY_Jets_ext1, DY_Jets_ext2, W_Jets.
Production timestamps: 260930_022240, 260930_022305, 260930_022328,
260930_022353 respectively. Input lists contain all available ROOT files.
The initial diagnostic takes 100 evenly spaced file indices per production;
exact processed paths are saved in out/prompt_heep_100files/*_files.txt.
DY components are combined with unit electron weights. No luminosity,
cross-section, pileup, or signed generator weights are applied. Results are
unweighted diagnostic efficiencies in this selected, sampled ntuple population.

Population:
  Reco raw pT>20 (threshold study) or >35 GeV (common-ID comparison).
  Nearest hard-process final-state gen electron with |eta_gen|<2.5
  within dR<0.001; one entry per reco electron, no opposite-charge pairing.
  EB |etaSC|<1.4442; EE 1.566<|etaSC|<2.5.
  No HEEP or isolation preselection. This is a prompt gen-matched sample,
  NOT a verified PF-isolated sample: no PF isolation or full particle/jet
  collection is available in these flat ntuples.
  Existing ntuple minLeptonsSkim=2 already biases the retained population,
  especially W+jets. These are not inclusive W-electron efficiencies.

Cut decisions are read from electron_heepV70Bitmap and electron_modHeepBitmap.
The stored pass flags are cross-checked against bitmap==4095.
Independent cut passing = passed that bit / selected electrons.
Cumulative passing = passed all bits through that bit / selected electrons.
N-minus-one failure = fails that bit but passes other 11 / passes other 11.
If the N-minus-one denominator is zero, the plot uses zero as a placeholder;
inspect cuts.csv denominator columns before interpreting an empty category.
Recovered/lost counts compare the two IDs on exactly the same electrons.
Full counters, cumulative counters, and per-cut transitions are stored in CSV.
SC acceptance and pT>35 are imposed in the main denominator, so their cuts
are already satisfied; pT>20 plots expose the nominal 35-GeV threshold.

Input plots use pT>35 and retain the sign of dEtaInSeed/dPerpIn/dPhiIn.
Four curves: all matched, full nominal pass, full modified pass,
and nominal bit2 fail / modified bit2 pass (not necessarily full-ID recovery).
Every curve is normalized to its own valid entries including flow bins.
Reject nonfinite, -999, and |value|>=1e30; dPerpIn can be unavailable.
flows.csv records valid entries and under/overflow. Raw ROOT histograms remain
unscaled. matching_same_subset compares absolute dEtaInSeed and dPerpIn on
the same electrons with both values valid; EE dPerpIn is not used by the ID.

No nominal isolation value is reconstructed from electron_dr03TkSumPt:
HEEP V70 uses dr03TkSumPtHEEP, which is not stored here. Likewise the branch
electron_modEcalHcalIso contains ECAL-only values and may lack matching product
references. Isolation cut effects therefore come from the authoritative VID
bits; modTrkIso plots display only the available modified input.
