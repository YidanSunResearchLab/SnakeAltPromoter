#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import subprocess
import json
import re
import glob
import sys
import csv
import shutil

def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Run a comprehensive pipeline for alternative promoter analysis.")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to the input directory containing paired-end FASTQ gz files (e.g., /path/to/fastqs/).")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output directory where results will be saved (e.g., /path/to/output/).")
    parser.add_argument("--organism", required=True, help="Reference genome assembly to use, created by the Genomesetup step. For example: 'hg38', 'dm6', 'ce3', or 'mm10'.")
    parser.add_argument("--genome_dir", required=True, help="Absolute path to the directory containing pre-generated genome files, created by the Genomesetup step (e.g., /absolute/path/to/genomes/).")
    parser.add_argument("--downsample_size", type=int, default=0, help="Number of valid pairs to downsample to for analysis. Set to 0 to disable downsampling (default: 0).")
    parser.add_argument("--fastqc",action="store_true",help="Enable this flag if needs fastqc.")
    parser.add_argument("--trim",action="store_true",help="Enable this flag if reads were trimmed using Trim Galore. If not set, the pipeline will use downsampled fastqs.")
    parser.add_argument("--threads", type=int, default=30, help="Number of CPU threads to use for parallel processing (default: 30).")
    parser.add_argument("--method", type=str, default="rnaseq", choices=["salmon", "proactiv", "dexseq", "cage", "rnaseq"], help="Which method to run: salmon / proactiv / dexseq / cage / rnaseq")
    parser.add_argument("--reads", type=str, default="paired", choices=["single", "paired"], help=" Reads are single-ended or paired: single / paired")
    parser.add_argument("--min_pFC", type=float, default=2.0, help="Additional threshold of minimum fold change of promoter activity for a promoter to be considered alternative promoter (default 2.0)")
    parser.add_argument("--max_gFC", type=float, default=1.5, help="Additional threshold of maximum fold change of gene expression for a promoter to be considered alternative promoter (default 1.5)")
    parser.add_argument("--lfcshrink", action="store_true", help="Enable log2 fold change shrinkage during differential analysis.")
    parser.add_argument("--sample_sheet", type=str, default=None, help="Path to sampleSheet.tsv file. Contains 'sampleName', 'condition', 'batch', and 'differential' columns."
                                                                        "If not provided, a default will be created automatically in the output directory named samplesheet.tsv.")
    parser.add_argument("--slurm", action="store_true", help="Use Snakemake native SLURM executor (--slurm).")
    parser.add_argument("--slurm-account", default=None, help="SLURM account; passed via --default-resources slurm_account=<...>.")
    parser.add_argument("--slurm-partition", default=None, help="SLURM partition; passed via --default-resources slurm_partition=<...>.")
    parser.add_argument("--set-resources", action="append", default=[], help="Per-rule resource override like '<rule>:slurm_partition=<PART>'. Repeatable. Mirrors Snakemake docs.")

    # Parse known arguments, capturing extra Snakemake args
    args, extra_args = parser.parse_known_args()

    # Assign variables from parsed arguments
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    organism = args.organism
    genome_dir = args.genome_dir
    downsample_size = args.downsample_size
    fastqc = args.fastqc
    trim = args.trim
    reads_choice = args.reads.lower()
    if reads_choice == "paired":
        reads = ["R1", "R2"]
    elif reads_choice == "single":
        reads = ["R1"]
    else:
        raise ValueError("Invalid --reads value. Must be 'single' or 'paired'.")



    threads = args.threads
    method = args.method
    min_pFC = args.min_pFC
    max_gFC = args.max_gFC
    lfcshrink = args.lfcshrink

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Determine script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Check if required files exist
    required_files = [
        f"{genome_dir}/organisms/{organism}/STARIndex",
        f"{genome_dir}/organisms/{organism}/Annotation.gtf",
        f"{genome_dir}/organisms/{organism}/Annotation/genes.bed",
        f"{genome_dir}/organisms/{organism}/Annotation/DEXSeq_flattened_exons.gff",
        f"{genome_dir}/organisms/{organism}/Annotation/proActiv_promoter_annotation.rds",
        f"{genome_dir}/organisms/{organism}/Annotation/genes_t2g.tsv"
    ]

    print("Checking required files...")  # Debug
    for f in required_files:
        print(f"Checking: {f}")  # Debug
        if not os.path.exists(os.path.abspath(f)):
            error_msg = (
                f"Dear user: Genome setup is not done. Required file '{f}' not found. "
                f"Since this is likely the first time running the pipeline for '{organism}', "
                "please run Genomesetup.py first with the path to the organism's FASTA file.\n"
            )
            sys.stderr.write(error_msg)
            sys.stderr.flush()  # Ensure message is written
            print("Exiting due to missing file.")  # Debug to stdout
            sys.exit(1)
    print("All required files found.")  # Debug


    def detect_samples(input_dir, reads_choice):
        """
        Return dict(sample -> {"R1": path, "R2": path (optional)})
        - For paired: require SAMPLE_R1.fastq.gz and/or SAMPLE_R2.fastq.gz
        - For single: accept SAMPLE_R1.fastq.gz OR SAMPLE.fastq.gz
        """
        samples_detected = {}

        if reads_choice == "paired":
            pat = re.compile(r"^(?P<sample>.+)_R(?P<read>[12])\.fastq\.gz$")
            for fn in os.listdir(input_dir):
                m = pat.match(fn)
                if not m:
                    continue
                sample = m.group("sample")
                read   = m.group("read")
                key    = f"R{read}"
                samples_detected.setdefault(sample, {})[key] = os.path.join(input_dir, fn)

            for sample, files in samples_detected.items():
                if "R1" not in files:
                    raise ValueError(f"Sample {sample} missing R1 file.")
                if "R2" not in files:
                    print(f"Warning: Sample {sample} missing R2 – treated as single-end.", file=sys.stderr)

        else:  # reads_choice == "single"
            pat_r1 = re.compile(r"^(?P<sample>.+)_R1\.fastq\.gz$")
            pat_plain = re.compile(r"^(?P<sample>.+)\.fastq\.gz$")
            for fn in os.listdir(input_dir):
                m1 = pat_r1.match(fn)
                m2 = pat_plain.match(fn) if not m1 else None
                if not (m1 or m2):
                    continue
                sample = (m1 or m2).group("sample")
                # If both SAMPLE_R1.fastq.gz and SAMPLE.fastq.gz exist, prefer the explicit R1
                cur = samples_detected.setdefault(sample, {})
                if m1:
                    cur["R1"] = os.path.join(input_dir, fn)
                elif "R1" not in cur:
                    cur["R1"] = os.path.join(input_dir, fn)

            if not samples_detected:
                raise ValueError(f"No single-end FASTQs found in {input_dir}. "
                                f"Expect files named SAMPLE.fastq.gz or SAMPLE_R1.fastq.gz.")

        return samples_detected



    def load_sample_sheet(sample_sheet_path):
        """Return (ordered list of samples, dict(sample->condition), dict(sample->batch), dict(sample->role), test_condition: str, control_condition: str)"""
        if not os.path.exists(sample_sheet_path):
            raise FileNotFoundError(f"sampleSheet.tsv not found: {sample_sheet_path}")

        sample_order = []
        cond_map = {}
        batch_map = {}
        role_map = {}

        with open(sample_sheet_path, newline="") as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            required = {"sampleName", "condition", "differential"}
            missing_cols = required - set(reader.fieldnames or [])
            if missing_cols:
                raise ValueError(f"sampleSheet.tsv is missing columns: {missing_cols}")

            has_batch = "batch" in reader.fieldnames
            role2conds = {"test": set(), "control": set()}


            for row in reader:
                sample = row["sampleName"].strip()
                cond   = row["condition"].strip()
                role   = row["differential"].strip().lower()
                if role not in {"test", "control"}:
                    raise ValueError(
                        f"Row for sample '{sample}' has invalid differential='{row['differential']}'. "
                        f"Allowed: 'test' or 'control'."
                    )
                sample_order.append(sample)
                cond_map[sample] = cond
                role_map[sample] = role
                role2conds[role].add(cond)

                # if no 'batch' column, default to 'batch1'
                if has_batch:
                    batch_map[sample] = (row["batch"].strip() or "batch1")
                else:
                    batch_map[sample] = "batch1"
        # Handle differential roles
        if role2conds["test"] and role2conds["control"]:
            # Standard two-condition case
            if len(role2conds["test"]) != 1 or len(role2conds["control"]) != 1:
                raise ValueError(
                    f"'differential' roles must each map to exactly ONE condition. "
                    f"Found: test -> {sorted(role2conds['test'])}, "
                    f"control -> {sorted(role2conds['control'])}"
                )
            test_condition    = next(iter(role2conds["test"]))
            control_condition = next(iter(role2conds["control"]))
        else:
            # Single-condition fallback
            only_cond = next(iter(role2conds["test"] or role2conds["control"]))
            test_condition = only_cond
            control_condition = only_cond
            print(f"[INFO] Only one condition detected ('{only_cond}'). Differential testing will be skipped.", file=sys.stderr)

        return sample_order, cond_map, batch_map, role_map, test_condition, control_condition

    # Main usage
    #sample_sheet_path = args.sample_sheet or os.path.join(input_dir, "sampleSheet.tsv")
    samples_detected = detect_samples(input_dir, reads_choice)

    # Determine sample sheet path
    if args.sample_sheet:
        sample_sheet_path = args.sample_sheet
    else:
        sample_sheet_path = os.path.join(output_dir, "samplesheet.tsv")

    # Auto-generate sampleSheet if missing
    if not os.path.exists(sample_sheet_path):
        os.makedirs(os.path.dirname(sample_sheet_path), exist_ok=True)
        with open(sample_sheet_path, "w", newline="") as tsvout:
            writer = csv.writer(tsvout, delimiter="\t")
            writer.writerow(["sampleName", "condition", "batch", "differential"])
            for sample in sorted(samples_detected.keys()):
                writer.writerow([sample, "wt", "batch1", "control"])
        print(f"[INFO] Created default sample sheet at: {sample_sheet_path}", file=sys.stderr)
        print("[INFO] Defaults used: condition=wt, batch=batch1, differential=control", file=sys.stderr)

    sample_order, condition_map, batch_map, role_map, test_condition, control_condition = load_sample_sheet(sample_sheet_path)


    # Check for missing samples
    missing = [s for s in sample_order if s not in samples_detected]
    if missing:
        raise ValueError(f"Samples in sampleSheet but not in input_dir: {missing}")

    # Build ordered samples_dict with condition and batch attached
    samples_dict = {
        s: {
            "R1": samples_detected[s]["R1"],
            "R2": samples_detected[s].get("R2"),
            "condition": condition_map[s],
            "batch": batch_map[s],
        }
        for s in sample_order
    }

    # Strings to pass to R, strictly ordered by sample_order
    conditions_str = ",".join(condition_map[s] for s in sample_order)
    batch_str      = ",".join(batch_map[s] for s in sample_order)

    # Print for verification
    print("Samples detected from directory & sampleSheet:")
    for sample, files in samples_dict.items():
        print(
            f"  {sample}: "
            f"R1={files['R1']}, "
            f"R2={files.get('R2', 'NA')}, "
            f"condition={files['condition']}, "
            f"batch={files['batch']}"
        )



    # Build Snakemake command
    snakemake_cmd = [
        "snakemake",
        "--snakefile", os.path.join(script_dir, "../rules/Snakealtpromoter.Snakefile"),
        "--printshellcmds",
        "--directory", output_dir,
        "--use-conda",
        "--conda-prefix", os.path.join(genome_dir, ".snakemake_conda"),
        "--config",
        f"input_dir={input_dir}",
        f"output_dir={output_dir}",
        f"organism={organism}",
        f"genome_dir={genome_dir}",
        f"downsample_size={downsample_size}",
        f"fastqc={fastqc}",
        f"trim={trim}",
        #f"samples={json.dumps(samples)}",
        f"reads={json.dumps(reads)}",
        f"threads={threads}",
        f"samples_dict={json.dumps(samples_dict)}",
        f"method={method}",
        #f"reads={','.join(reads)}",
        f"test_condition={test_condition}",
        f"control_condition={control_condition}",
        f"max_gFC={max_gFC}",
        f"min_pFC={min_pFC}",
        f"lfcshrink={lfcshrink}",
        #f"conditions={conditions_str}",
        f"batch={batch_str}",
        "--cores", str(threads),
    ]

    # SLURM options
    if args.slurm:
        # ensure slurm client available
        if shutil.which("sbatch") is None:
            print("[ERROR] --slurm requested but 'sbatch' not found. Run on a SLURM node.", file=sys.stderr)
            sys.exit(2)

        # enable native slurm executor
        snakemake_cmd.append("--slurm")

        # add default resources (account, partition)
        default_resources = []
        if args.slurm_account:
            default_resources.append(f"slurm_account={args.slurm_account}")
        if args.slurm_partition:
            default_resources.append(f"slurm_partition={args.slurm_partition}")

        if default_resources:
            snakemake_cmd += ["--default-resources"] + default_resources

        # per-rule overrides (e.g. <rule>:slurm_partition=highmem)
        for res in args.set_resources:
            snakemake_cmd += ["--set-resources", res]

    # Append any extra Snakemake arguments
    if extra_args:
        snakemake_cmd.extend(extra_args)

    # Generate unique log file name
    log_pattern = os.path.join(output_dir, "snakealtpromoter_run_*.log")
    log_number = max([int(re.search(r"snakealtpromoter_run_(\d+)\.log", f).group(1))
                    for f in glob.glob(log_pattern)], default=0) + 1
    log_file_path = os.path.join(output_dir, f"snakealtpromoter_run_{log_number}.log")

    # Execute Snakemake with logging
    with open(log_file_path, "w") as log_file:
        process = subprocess.Popen(
            snakemake_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        for line in process.stdout:
            print(line, end="")
            log_file.write(line)
        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, snakemake_cmd)


if __name__ == "__main__":
    main()
