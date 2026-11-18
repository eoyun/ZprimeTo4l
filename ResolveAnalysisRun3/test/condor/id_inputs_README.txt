Electron ID input distributions for HToAATo4L signal

Run from this directory:
  source env.sh
  python3 -u plot_id_inputs.py
Redraw without reading CRAB ntuples:
  python3 plot_id_inputs.py --plots-only

Inputs: existing filelists/HToAATo4L_H*A*_22.txt, all listed files.
Tree: resolvedNtuplizer/Events.
Outputs: out/id_inputs/<sample>.root (unweighted raw histograms),
         out/id_inputs/counts.csv (selected/valid/invalid/underflow/overflow),
         imgs/id_inputs/*.png and *.pdf.

Selection follows the existing HeepEffiCore.h denominator:
  Greedy nearest opposite-charge gen electron pairing (not ancestry-based).
  Both gen electrons |eta|<2.5 and pair leading gen pT>20 GeV.
  Reco electron raw pT>10 GeV and nearest selected gen electron dR<0.001.
  Each reco electron counted once. No HEEP or modified-HEEP requirement.
  EB: |etaSC|<1.4442. EE: 1.566<|etaSC|<2.5.
  Unit weights; no luminosity/cross-section/generator-weight scaling.

Branches: electron_dPerpIn, electron_dEtaInSeed, electron_sigmaIeta.
Matching variables retain their sign (not absolute values).
Nonfinite values, -999 and |value|>=1e30 are excluded and counted invalid.
In particular, dPerpIn without a usable extra track is not a measurement.
The dPerpIn valid subset may differ from the other variables' subsets.

Each plotted histogram is normalized to all its valid entries, including
underflow and overflow. Out-of-range entries are not folded into edge bins;
their counts are in counts.csv and the raw ROOT histogram flow bins.
Matching-variable displayed ranges: [-0.03, 0.03], 120 bins.
Sigma displayed ranges: EB [0, 0.03], EE [0, 0.10], 120 bins.
Plots use log-y. Each sample has a 3-column (variable), 2-row (EB/EE) page.
compare_H*: A-mass overlays at fixed H mass, EB and EE side by side.
compare_A*: H-mass overlays at fixed A mass, EB and EE side by side.
