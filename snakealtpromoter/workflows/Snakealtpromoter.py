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
import yaml

def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def detect_samples_bam(input_dir):
    """
    BAM-based detection for precomputed alignments.

    Expected directory layout:

        <input_dir>/<sample>/genome_alignments_out/
            genome_accepted_hits.bam
            SJ.out.tab

    Returns:
        {
            sample_name: {
                "bam": <path_to_genome_accepted_hits.bam>,
                "sj":  <path_to_SJ.out.tab>,
            },
            ...
        }
    """
    samples_detected = {}

    for sample in os.listdir(input_dir):
        sample_dir = os.path.join(input_dir, sample)
        if not os.path.isdir(sample_dir):
            continue

        ga_dir = os.path.join(sample_dir, "genome_alignments_out")
        if not os.path.isdir(ga_dir):
            continue

        bam_path = os.path.join(ga_dir, "genome_accepted_hits.bam")
        bai_path = os.path.join(ga_dir, "genome_accepted_hits.bam.bai")
        sj_path  = os.path.join(ga_dir, "SJ.out.tab")

        if os.path.exists(bam_path) and os.path.exists(sj_path) and os.path.exists(bai_path):
            samples_detected[sample] = {
                "bam": bam_path,
                "sj":  sj_path,
                "bai": bai_path,
            }
        else:
            print(
                f"[WARN] sample '{sample}' is missing genome_accepted_hits.bam or SJ.out.tab or genome_accepted_hits.bam.bai; skipped.",
                file=sys.stderr
            )

    return samples_detected


def detect_samples_with_fallback(input_dir, reads_choice):
    """
    Wrapper around the original FASTQ-based detect_samples(..).

    Behavior:
      1) Try FASTQ-based detection first (original behavior).
         - If some FASTQ samples are found, use them directly.
         - If single-end mode raises the 'No single-end FASTQs found' error,
           treat it as 'no FASTQ samples found' and fall back.
      2) If no FASTQ samples are found, fall back to BAM-based detection
         with layout:
             <input_dir>/<sample>/genome_alignments_out/
                 genome_accepted_hits.bam
                 SJ.out.tab
    """
    # ---- Step 1: try original FASTQ detection ----
    try:
        fastq_samples = detect_samples(input_dir, reads_choice)
    except ValueError as e:
        msg = str(e)
        # Only swallow the specific "no single-end FASTQs" error;
        # any other ValueError should still be raised.
        if "No single-end FASTQs found" in msg:
            print(
                f"[INFO] FASTQ-based detection raised '{msg}'. "
                f"Will try BAM-based detection instead.",
                file=sys.stderr,
            )
            fastq_samples = {}
        else:
            # Different error → propagate
            raise

    if fastq_samples:
        print("[INFO] Detected FASTQ-based samples; using FASTQ mode.", file=sys.stderr)
        return fastq_samples

    # ---- Step 2: fall back to BAM layout ----
    print(
        "[INFO] No FASTQ samples found; trying BAM-based layout "
        "under <sample>/genome_alignments_out/ ...",
        file=sys.stderr,
    )
    bam_samples = detect_samples_bam(input_dir)

    if not bam_samples:
        raise ValueError(
            f"No samples found in {input_dir} in either FASTQ mode or BAM mode.\n"
            f"Expected either FASTQs (*.fastq.gz) or BAM layout "
            f"<sample>/genome_alignments_out/."
        )

    print("[INFO] Detected BAM-based samples; using BAM mode.", file=sys.stderr)
    return bam_samples

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
    parser.add_argument("--star_precomputed", help="Path to precomputed STAR genome index directory.")

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
    star_precomputed = args.star_precomputed
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

    def load_sample_sheet(sample_sheet_path):
        """Return (ordered list of samples, dict(sample->condition), dict(sample->batch), dict(sample->role), test_condition: str, control_condition: str)"""
        
        if not os.path.exists(sample_sheet_path):
            raise FileNotFoundError(f"sampleSheet.tsv not found: {sample_sheet_path}")
        
        sheet_name = os.path.splitext(os.path.basename(sample_sheet_path))[0]

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

        return sample_order, cond_map, batch_map, role_map, test_condition, control_condition, sheet_name

    # Main usage
    # sample_sheet_path = args.sample_sheet or os.path.join(input_dir, "sampleSheet.tsv")

    # Use FASTQ-based detection first; if none found, fall back to BAM layout
    samples_detected = detect_samples_with_fallback(input_dir, reads_choice)
    EXCLUDE_SAMPLES = {"RISK_277_S88", "1075-DLPFC", "2007-DLPFC", "41_120416"}

    for x in EXCLUDE_SAMPLES:
        if x in samples_detected:
            samples_detected.pop(x, None)
            print(f"[INFO] Excluding sample from detection: {x}", file=sys.stderr)



    # Determine sample sheet path
    if args.sample_sheet:
        sample_sheet_path = args.sample_sheet
    else:
        sample_sheet_path = os.path.join(output_dir, "samplesheet.tsv")

    # Auto-generate a TEMPLATE sampleSheet if missing (no grouping inferred)
    if not os.path.exists(sample_sheet_path):
        os.makedirs(os.path.dirname(sample_sheet_path), exist_ok=True)
        detected_sorted = sorted(samples_detected.keys())

        with open(sample_sheet_path, "w", newline="") as tsvout:
            writer = csv.writer(tsvout, delimiter="\t")
            writer.writerow(["sampleName", "condition", "batch", "differential"])
            for sample in detected_sorted:
                # Template defaults: no differential plan
                writer.writerow([sample, "overall", "batch1", "control"])

        print(f"[INFO] Created TEMPLATE sample sheet at: {sample_sheet_path}", file=sys.stderr)
        print("[INFO] This template disables differential by default.", file=sys.stderr)
        print("[INFO] Provide your own sampleSheet (via --sample_sheet) with test/control + condition labels to run differential.", file=sys.stderr)

    sample_order, condition_map, batch_map, role_map, test_condition, control_condition, sheet_name = load_sample_sheet(sample_sheet_path)

    # ---- Validate sampleSheet entries (only those listed) exist in detected samples ----
    missing = [s for s in sample_order if s not in samples_detected]
    if missing:
        raise ValueError(f"Samples in sampleSheet but not in input_dir: {missing}")

    # ---- Define ALL samples from detection (pipeline runs these) ----
    all_samples = sorted(samples_detected.keys())

    # ---- Define DIFF samples from sampleSheet (only these used for differential) ----
    diff_samples = [
        s for s in sample_order
        if (s in samples_detected) and (role_map.get(s) in {"test", "control"})
    ]

    # Build samples_dict for ALL samples; condition/batch come from sampleSheet if present, else defaults
    samples_dict = {}
    for s in all_samples:
        sd = samples_detected[s]

        entry = {
            "condition": condition_map.get(s, "overall"),
            "batch": batch_map.get(s, "batch1"),
        }

        # FASTQ mode
        if "R1" in sd:
            entry["R1"] = sd["R1"]
        if "R2" in sd:
            entry["R2"] = sd["R2"]

        # BAM mode (expected from detect_samples_bam)
        if "bam" in sd:
            entry["bam"] = sd["bam"]
        if "sj" in sd:
            entry["sj"] = sd["sj"]
        if "bai" in sd:
            entry["bai"] = sd["bai"]

        samples_dict[s] = entry

    # Strings to pass to R for ALL samples (if you still need them)
    all_conditions_str = ",".join(samples_dict[s]["condition"] for s in all_samples)
    all_batch_str      = ",".join(samples_dict[s]["batch"] for s in all_samples)

    # Strings to pass to R for DIFF samples (what differential should use)
    diff_conditions_str = ",".join(samples_dict[s]["condition"] for s in diff_samples)
    diff_batch_str      = ",".join(samples_dict[s]["batch"] for s in diff_samples)
    diff_samples_str    = ",".join(diff_samples)

    print("Samples detected (ALL):")
    for s in all_samples:
        files = samples_dict[s]
        print(
            f"  {s}: "
            f"R1={files.get('R1','NA')}, "
            f"R2={files.get('R2','NA')}, "
            f"bam={files.get('bam','NA')}, "
            f"sj={files.get('sj','NA')}, "
            f"bai={files.get('bai','NA')}, "
            f"condition={files['condition']}, "
            f"batch={files['batch']}"
        )

    print("Differential samples (DIFF subset):", diff_samples, file=sys.stderr)

    # ---- Build CONFIG dict ----
    CONFIG = {
        "input_dir": input_dir,
        "output_dir": output_dir,
        "organism": organism,
        "genome_dir": genome_dir,
        "downsample_size": downsample_size,
        "fastqc": fastqc,
        "trim": trim,
        "reads": reads,
        "threads": threads,
        "samples_dict": samples_dict,
        "method": method,
        "test_condition": test_condition,
        "control_condition": control_condition,
        "sheet_name": sheet_name,
        "max_gFC": max_gFC,
        "min_pFC": min_pFC,
        "lfcshrink": lfcshrink,
        "all_batch_str": all_batch_str,
        "diff_batch_str": diff_batch_str,
        "diff_samples_str": diff_samples_str,
        "diff_conditions_str": diff_conditions_str,
        "star_precomputed": args.star_precomputed,
    }

    # ---- Write CONFIG to YAML ----
    configfile_path = os.path.join(output_dir, "Snakealtpromoter_config.yaml")
    os.makedirs(output_dir, exist_ok=True)
    with open(configfile_path, "w") as fh:
        yaml.safe_dump(CONFIG, fh, sort_keys=False)


    snakemake_cmd = [
        "snakemake",
        "--snakefile", os.path.join(script_dir, "../rules/Snakealtpromoter.Snakefile"),
        "--printshellcmds",
        "--directory", output_dir,
        "--use-conda",
        "--conda-prefix", os.path.join(genome_dir, ".snakemake_conda"),
        "--configfile", configfile_path,
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
