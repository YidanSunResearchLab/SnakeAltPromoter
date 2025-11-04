#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import subprocess
import glob
import re
import sys

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Execute the AltPromoterFlow genome setup pipeline to generate genome indices and annotations.")
    parser.add_argument("-o", "--output_dir", required=True, help="Absolute path to the output directory for storing generated genome files.")
    parser.add_argument("--organism", required=True, help="Reference genome assembly (e.g., 'hg38', 'mm10', 'dm6').")
    parser.add_argument("--organism_fasta", required=True, help="Path to the organism FASTA file with 'chr' prefix (e.g., /path/to/hg38.fa).")
    parser.add_argument("--genes_gtf", required=True, help="Path to the GTF file for gene annotations (e.g., /path/to/gencode.v38.annotation.gtf).")
    parser.add_argument("--threads", type=int, default=16, help="Number of CPU threads for parallel processing (default: 16).")

    # Parse known arguments, capturing extra Snakemake args
    args, extra_args = parser.parse_known_args()

    # Ensure output directory exists
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Determine script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Validate input files
    if not os.path.exists(args.organism_fasta):
        sys.stderr.write(f"Error: FASTA file '{args.organism_fasta}' does not exist.\n")
        sys.exit(1)
    if not os.path.exists(args.genes_gtf):
        sys.stderr.write(f"Error: GTF file '{args.genes_gtf}' does not exist.\n")
        sys.exit(1)

    # Build Snakemake command
    snakemake_cmd = [
        "snakemake",
        #"--snakefile", os.path.join(script_dir, "genomesetup.smk"),
        "--snakefile", os.path.abspath(os.path.join(script_dir, "../rules/Genomesetup.Snakefile")),
        "--printshellcmds",
        "--directory", output_dir,
        "--use-conda",
        "--conda-prefix", os.path.join(output_dir, ".snakemake_conda"),
        "--config",
        f"output_dir={output_dir}",
        f"organism={args.organism}",
        f"organism_fasta={args.organism_fasta}",
        f"genes_gtf={args.genes_gtf}",
        f"threads={args.threads}",
        "--cores", str(args.threads),
    ]

    # Append any extra Snakemake arguments
    if extra_args:
        snakemake_cmd.extend(extra_args)

    # Generate unique log file name
    log_pattern = os.path.join(output_dir, f"genomesetup_run_{args.organism}_*.log")
    log_number = max([int(re.search(r"genomesetup_run_.+_(\d+)\.log", os.path.basename(f)).group(1))
                     for f in glob.glob(log_pattern)], default=0) + 1
    log_file_path = os.path.join(output_dir, f"genomesetup_run_{args.organism}_{log_number}.log")

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