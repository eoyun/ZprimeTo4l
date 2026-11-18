# CRAB output merging

Each sample name is preserved, including `ext`, `ext1`, and `ext2`. Data eras
are `2022C`/`2022D`; the current MC campaign is `2022preEE`. Existing CRAB
submission scripts are not modified.

## Initial execution (2026-10-08)

- Smoke DAG `work/20261008_224528_522399`: cluster `9733694`, succeeded.
  ZZ's two files merged with ROOT 6.26.11 and verified tree entries 36/1.
- Production DAG `work/20261008_224630_053698`: cluster `9733700`, submitted.
  27 data/background samples, 48,411 input files, 193 nodes, at most four jobs.
  This is an active run; do not submit the same DAG again.
- Signal status refresh: `work/20261008_225531_969170/status` contains current
  CRAB responses. Signals with failed jobs are excluded pending CRAB recovery.

From the original submission folder:

```bash
python3 merge/status_merge.py merge/work/20261008_224630_053698
condor_q 9733700
```

## Prepare and inspect

```bash
python3 merge/prepare_merge.py
```

Inspect the printed run directory's `manifest.json`. Only completed samples
with the expected number of files enter the DAG. Saved logs are historical;
to query current status, initialize the CRAB client and use:

```bash
merge/refresh_status.sh --samples HToAATo4L_H250A1_22
```

`--samples` accepts exact names. Unrecognized eras require `era_overrides`
in `config.json`. Multiple submission timestamps are rejected. Inputs come
only from numbered subdirectories of the chosen task; logs and other ROOT
outputs are excluded. No original and extension datasets are combined.

## Smoke test and submit

```bash
python3 merge/prepare_merge.py --smoke --samples ZZ
# Submit the DAG path printed by the command:
condor_submit_dag -maxjobs 4 /absolute/path/to/work/run/merge.dag
```

A smoke run uses at most two files and writes under its own run directory.
Submit a normal prepared DAG with the same command. DAGMan retries failed
nodes twice and starts final merging only after all partial merges succeed.
The default maximum of four jobs limits storage load.

Inspect progress without contacting the scheduler:

```bash
python3 merge/status_merge.py /absolute/path/to/work/run
```

Workers require `/pnfs`, `/u`, CVMFS, and Apptainer/Singularity. ROOT is loaded
from CVMFS CMSSW_13_0_13 in the CMS EL8 container. Runtime scripts are transferred
by Condor; manifests, lists, intermediate outputs and results use shared storage.
The default disk request is 50 GiB per job; adjust after inspecting file sizes.

## Verification and recovery

Each node checks object and tree branch schemas, total tree entries, and named
normalization/cutflow histograms before publishing. `hadd` merges both trees
and histograms and does not skip corrupt inputs. A JSON sidecar records inputs,
ROOT version, tree entries, and normalization values. Input-list hashes are in
the manifest. Histogram bin errors and every analysis histogram are not compared
individually. A test on real inputs is required before large production runs.

Results are in `merge/results/<era>/<sample>.root`; intermediate files remain
under the run directory. Outputs are never overwritten. Publication uses a
destination-side staging file and an atomic hard link. On filesystems without
hard-link support this fails rather than silently changing publication behavior.

Read `logs/<node>.out`, `logs/<node>.err`, `logs/condor.log`, and DAGMan logs
for failures. Resume the existing DAG using DAGMan's rescue DAG; preparing a new
run creates new intermediate paths. A worker killed abruptly may leave an
`*.lock`; inspect the job and destination before removing that reservation.
