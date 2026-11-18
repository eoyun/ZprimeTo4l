DY/W cumulative cutflow only

source env.sh
python3 plot_dyw_cutflow.py
python3 plot_dyw_cutflow.py --plots-only

Outputs: imgs/dyw_cutflow_only/{DY,W}_cutflow_pt{20,35}.{png,pdf}
Each figure: prompt matched (top), unmatched (bottom); EB left, EE right.
Only cumulative cutflows are drawn. Blue dashed: nominal HEEP;
red solid: modified HEEP. Baseline is normalized to 1 for each
population/region. Legend gives final passing count / baseline count.
Raw counts, fraction and loss at each step: out/dyw_cutflow_only/cutflow.csv.

Both populations use the SAME gen matching definition: nearest stored
hard-process final-state electron with |gen eta|<2.5, matched if dR<0.1.
No eligible gen electron means unmatched. This replaces the previous
prompt-only diagnostic's dR<0.001 definition. W/Z ancestry is not available;
prompt matched is shorthand for matched to the stored hard-process electrons.
Unmatched is a fake-candidate proxy, not proof of hadronic jet origin.

Fiducial: raw pT>20 or >35 GeV; EB |etaSC|<1.4442;
EE 1.566<|etaSC|<2.5. Existing two-lepton ntuple skim remains.
No ID/isolation preselection. Unit electron weights; DY extensions pooled.
Exact same 400 file manifests as out/prompt_heep_100files are reused and copied.

Use pT>35 figures to compare the modified cuts at a common pT threshold.
pT>20 figures include the nominal 35 vs modified 20 GeV ID threshold effect.
The pT/SC eta cuts are already satisfied by the >35 fiducial baseline.
Each stage includes all earlier cuts; effects depend on that ordering.
The script checks VID bits against stored pass flags, monotonicity, baseline
equality, and every unmatched cumulative count against the prior diagnostic.
