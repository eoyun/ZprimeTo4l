DY/W gen-unmatched reconstructed-electron ID passing rate

Run:
  source env.sh
  python3 -u compare_prompt_heep.py --population unmatched --match-dr 0.1 --max-files 100
Redraw: add --plots-only. Use --max-files 0 to process all listed files.
Outputs: out/heep_unmatched_dr0p1_100files and imgs/heep_unmatched_dr0p1_100files.

Same 100 evenly spaced files per production as the previous prompt comparison:
DY_Jets, DY_Jets_ext1, DY_Jets_ext2, W_Jets. DY components pooled with unit
electron weights. No luminosity/cross-section or signed MC weighting.

Gen reference: stored hard-process final-state electrons, |gen eta|<2.5.
Mother PDG ID / full ancestry is NOT stored. In these DY/W samples this is
used as a proxy for W/Z electrons, NOT a verified W/Z ancestry selection.
For each reco electron, find the closest eligible gen electron in eta/phi.
Matched means dR<0.1; keep dR>=0.1. If no eligible gen electron exists,
the reco candidate is unmatched and is included. No charge or gen-pT cut.

Reco baseline: raw pT>20 GeV, |etaSC|<1.4442 (EB), or
1.566<|etaSC|<2.5 (EE). Also report pT>35 GeV to compare IDs on common phase
space. No HEEP/loose-ID/isolation requirement before forming denominator.
Existing ntuple minLeptonsSkim=2 remains in effect.

Rate = unmatched reco electron candidates passing ID / unmatched reco
electron candidates in the same baseline. This is a reconstructed-electron
candidate fake-rate proxy, NOT a per-jet misidentification probability: jets
are not the denominator and no reco-jet association is available. Includes
nonprompt/conversion and any residual prompt matching failures; cannot prove
every unmatched electron is a hadronic jet fake from this flat ntuple.

Stored full-ID booleans are checked against the corresponding 12-bit maps.
summary.csv and rates.csv contain denominators, pass counts, nominal-only /
modified-only transitions and unweighted passing rates. rates.csv includes
68.27% Clopper-Pearson intervals, which describe counting uncertainty only.
cuts.csv: individual, cumulative, N-minus-one counts and cut transitions.
Input distributions: pT>35; each curve normalized to its own valid entries.
flows.csv records histogram flow counts. See prompt_heep_README.txt for the
shared cut-plot and missing-value conventions.
DY/W_passing_rate_vs_pt plots use variable pT bins and binomial intervals.
Exact pT=35 is included in the [35,40) differential bin, while the integrated
pT>35 summary strictly excludes that boundary. Above 1000 GeV is kept in
ROOT histogram overflow and in the integrated rate, not the plotted pT range.
