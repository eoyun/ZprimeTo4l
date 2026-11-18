#!/usr/bin/env python3
"""Freeze CRAB outputs into independent sample/era merge DAGs. Never submit implicitly."""
import argparse
from collections import Counter
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def assignment(text, key):
    match = re.search(re.escape(key) + r"\s*=\s*['\"]([^'\"]+)['\"]", text)
    if not match:
        raise ValueError(f"Missing setting: {key}")
    return match[1]


def era_for(sample, dataset, overrides):
    if sample in overrides:
        return overrides[sample]
    match = re.search(r"/Run(20\d\d[A-Z])-", dataset)
    if match:
        return match[1]
    if "Run3Summer22MiniAOD" in dataset or "Run3Summer22_MiniAOD" in dataset:
        return "2022preEE"
    raise ValueError(f"Unknown era for {sample}; set era_overrides")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.json")
    parser.add_argument("--samples", nargs="+", help="Exact sample names; ext remains separate")
    parser.add_argument("--refresh-status", action="store_true", help="Query CRAB for current status")
    parser.add_argument("--smoke", action="store_true", help="Separate test output using first two files")
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    if config["chunk_size"] < 1 or config["max_jobs"] < 1:
        raise ValueError("chunk_size and max_jobs must be positive")
    run = HERE / "work" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    for child in ["lists", "jobs", "logs", "partials", "status"]:
        (run / child).mkdir(parents=True, exist_ok=True)
    output_root = Path(config["output_root"])
    if not output_root.is_absolute():
        output_root = HERE / output_root
    if args.smoke:
        output_root = run / "smoke_results"
    rows = [line.split() for filename in ["inputdata.dat", "inputdataMC.dat"]
            for line in (BASE / filename).read_text().splitlines() if line.strip()]
    names = {row[0] for row in rows}
    if args.samples and not set(args.samples) <= names:
        raise ValueError(f"Unknown samples: {set(args.samples) - names}")
    dag, summary = [], []

    def add_job(node, paths, output):
        input_list = run / "lists" / (node + ".txt")
        input_list.write_text("\n".join(map(str, paths)) + "\n")
        job = {"input_list": str(input_list), "output": str(output),
               "cmssw_release": config["cmssw_release"], "scram_arch": config["scram_arch"]}
        job_path = run / "jobs" / (node + ".json")
        job_path.write_text(json.dumps(job, indent=2) + "\n")
        dag.append(f'JOB {node} {HERE / "hadd.sub"}')
        variables = {"runner": str(HERE / "run_hadd.sh"), "checker": str(HERE / "check_root.py"),
                     "setup": str(HERE / "setup_root.sh"), "job": str(job_path),
                     "run_dir": str(run), "node": node,
                     "memory_mb": config["request_memory_mb"], "disk_mb": config["request_disk_mb"]}
        # These paths are generated here, not shell command fragments.
        if any('"' in str(v) or "\n" in str(v) for v in variables.values()):
            raise ValueError("Unsupported quote/newline in DAG variable")
        dag.append("VARS " + node + " " + " ".join(f'{k}="{v}"' for k, v in variables.items()))
        dag.append(f'RETRY {node} {config["retries"]}')

    for sample, dataset in rows:
        if args.samples and sample not in args.samples:
            continue
        record = {"sample": sample, "dataset": dataset}
        summary.append(record)
        try:
            era = era_for(sample, dataset, config["era_overrides"])
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", era):
                raise ValueError("Invalid era path component")
            record["era"] = era
            directory = BASE / "crab" / sample
            project = directory / "crab_projects" / ("crab_resolvedMuCR_" + sample)
            settings = next(p for p in [directory / "submit_crab.py", directory / "submitMC_crab.py",
                                       directory / "submitSig_crab.py"] if p.exists()).read_text()
            if assignment(settings, "config.Data.inputDataset") != dataset:
                raise ValueError("Sample list and submitted dataset disagree")
            tag = assignment(settings, "config.Data.outputDatasetTag")
            log = (project / "crab.log").read_text(errors="replace")
            timestamps = list(dict.fromkeys(re.findall(r"Task name:\s+(\d{6}_\d{6}):", log)))
            if len(timestamps) != 1:
                raise ValueError(f"Ambiguous submission timestamps: {timestamps}")
            timestamp = timestamps[0]
            source = Path(config["storage_root"]) / dataset.split("/")[1] / tag / timestamp
            record.update(timestamp=timestamp, source=str(source))
            counts = re.findall(r"Jobs status:\s+finished\s+100\.0%\s+\((\d+)/(\d+)\)", log)
            latest_job_status = re.findall(r"Jobs status:\s+([^\n]+)", log)
            completed = bool(counts and latest_job_status and "100.0%" in latest_job_status[-1]
                             and "finished" in latest_job_status[-1])
            expected = int(counts[-1][1]) if completed else None
            record["status_source"] = "saved CRAB log"
            if args.refresh_status:
                result = subprocess.run(["crab", "status", "--json", "-d", str(project)],
                                        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        timeout=180)
                (run / "status" / (sample + ".txt")).write_text(result.stdout)
                if result.returncode:
                    raise ValueError(f"CRAB status failed; see status/{sample}.txt")
                # CRAB status --json emits a job-id -> job-detail JSON mapping among text lines.
                mappings = []
                decoder = json.JSONDecoder()
                for match in re.finditer(r"\{", result.stdout):
                    try:
                        value, _ = decoder.raw_decode(result.stdout[match.start():])
                    except ValueError:
                        continue
                    if isinstance(value, dict) and value and all(
                            str(k).isdigit() and isinstance(v, dict) and "State" in v
                            for k, v in value.items()):
                        mappings.append(value)
                if not mappings:
                    raise ValueError("Cannot parse current CRAB job states; refusing merge")
                states = max(mappings, key=len)
                record["job_states"] = dict(Counter(v["State"] for v in states.values()))
                completed = all(v["State"] == "finished" for v in states.values())
                expected = len(states) if completed else None
                record["status_source"] = "live CRAB status"
            record["expected_jobs"] = expected
            files = []
            for subdir in sorted(source.iterdir()):
                if subdir.name.isdigit() and subdir.is_dir():
                    files.extend(p for p in subdir.iterdir() if re.fullmatch(r"hists_\d+\.root", p.name))
            files.sort(key=lambda p: int(re.search(r"_(\d+)\.root$", p.name)[1]))
            ids = [int(re.search(r"_(\d+)\.root$", p.name)[1]) for p in files]
            record["input_count"] = len(files)
            if not files or len(ids) != len(set(ids)):
                raise ValueError("Empty input or duplicate CRAB job IDs")
            if not args.smoke and (not completed or len(files) != expected):
                raise ValueError("Completion unconfirmed or input count differs from finished jobs")
            if args.smoke:
                files = files[:2]
            all_inputs = run / "lists" / (sample + "_all.txt")
            all_inputs.write_text("\n".join(map(str, files)) + "\n")
            record["input_list"] = str(all_inputs)
            record["input_list_sha256"] = hashlib.sha256(all_inputs.read_bytes()).hexdigest()
            final_output = output_root / era / (sample + ".root")
            if final_output.exists():
                raise ValueError(f"Output already exists: {final_output}")
            record["output"] = str(final_output)
            size = config["chunk_size"]
            chunks = [files[i:i + size] for i in range(0, len(files), size)]
            if len(chunks) == 1:
                add_job(sample + "_final", files, final_output)
            else:
                parents, partials = [], []
                for i, chunk in enumerate(chunks):
                    node = f"{sample}_part_{i:04d}"
                    partial = run / "partials" / sample / f"part_{i:04d}.root"
                    add_job(node, chunk, partial)
                    parents.append(node)
                    partials.append(partial)
                final = sample + "_final"
                add_job(final, partials, final_output)
                dag.append("PARENT " + " ".join(parents) + " CHILD " + final)
            record["state"] = "ready"
            record["merge_jobs"] = len(chunks) + (len(chunks) > 1)
        except (ValueError, OSError, StopIteration, subprocess.TimeoutExpired) as error:
            record.update(state="blocked", reason=str(error))
        print(sample, record["state"], record.get("input_count", "?"), record.get("reason", ""), flush=True)
    (run / "merge.dag").write_text("\n".join(dag) + "\n")
    (run / "manifest.json").write_text(json.dumps({"config": config, "smoke": args.smoke,
                                                  "samples": summary}, indent=2) + "\n")
    print("Run directory:", run)
    print("Ready:", sum(r["state"] == "ready" for r in summary), "/", len(summary))
    if dag:
        print(f"Submit: condor_submit_dag -maxjobs {config['max_jobs']} {run / 'merge.dag'}")


if __name__ == "__main__":
    main()
