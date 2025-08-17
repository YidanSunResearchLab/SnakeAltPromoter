import os
import argparse
import subprocess
import json
import re
import glob
import sys
import csv


def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Run a comprehensive pipeline integrating HiC-Pro|FitHiChIP, and Chromap|MAPS for HiChIP/PLAC-Seq data analysis.")
    parser.add_argument("-p", "--pipeline", choices=["Maps", "Fithichip", "Hichipper", "Hicdcplus", "All"], default="Maps", help="Software used to perform the data analysis")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to the input directory containing paired-end FASTQ gz files (e.g., /path/to/fastqs/).")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output directory where results will be saved (e.g., /path/to/output/).")
    parser.add_argument("--organism", required=True, help="Reference genome assembly to use, created by the Genomesetup step. For example: 'hg38', 'dm6', 'ce3', or 'mm10'.")
    parser.add_argument("--genome_dir", required=True, help="Absolute path to the directory containing pre-generated genome files, created by the Genomesetup step (e.g., /absolute/path/to/genomes/).")
    parser.add_argument("--downsample_size", type=int, default=0, help="Number of valid pairs to downsample to for analysis. Set to 0 to disable downsampling (default: 0).")
    parser.add_argument("--threads", type=int, default=30, help="Number of CPU threads to use for parallel processing (default: 30).")
    parser.add_argument("--test_condition", default="NA", help="Referece condition in differential analysis (e.g., --test_condition 'Disease'). Default: NA.")
    parser.add_argument("--control_condition", default="NA", help="Baseline condition in differential analysis (e.g., --control_condition 'Healthy'). Default: NA.")
    parser.add_argument("--method", type=str, default="rnaseq", choices=["salmon", "proactiv", "dexseq", "cage", "rnaseq"], help="Which method to run: salmon / proactiv / dexseq / cage / rnaseq")
    parser.add_argument("--reads", type=str, default="paired", choices=["single", "paired"], help=" Reads are single-ended or paired: single / paired")
    parser.add_argument("--min_pFC", type=float, default=2.0, help="Additional threshold of minimum fold change of promoter activity for a promoter to be considered alternative promoter (default 2.0)")
    parser.add_argument("--max_gFC", type=float, default=1.5, help="Additional threshold of maximum fold change of gene expression for a promoter to be considered alternative promoter (default 1.5)")
    # Parse known arguments, capturing extra Snakemake args
    args, extra_args = parser.parse_known_args()

    # Assign variables from parsed arguments
    pipeline = args.pipeline
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    organism = args.organism
    genome_dir = args.genome_dir
    downsample_size = args.downsample_size
    reads_choice = args.reads.lower()
    if reads_choice == "paired":
        reads = ["R1", "R2"]
    elif reads_choice == "single":
        reads = ["R1"]
    else:
        raise ValueError("Invalid --reads value. Must be 'single' or 'paired'.")

    

    threads = args.threads
    method = args.method
    test_condition = args.test_condition
    control_condition = args.control_condition
    min_pFC = args.min_pFC
    max_gFC = args.max_gFC

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

    def detect_samples(input_dir):
        """
        Return dict(sample -> {"R1": path, "R2": path (optional)})
        """
        pattern = re.compile(r"^(?P<sample>.+)_R(?P<read>[12])\.fastq\.gz$")
        samples_detected = {}

        for fn in os.listdir(input_dir):
            m = pattern.match(fn)
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

        return samples_detected


    def load_conditions(sample_sheet_path):
        """Return ordered list of sample names and dict(sampleName -> condition)"""
        if not os.path.exists(sample_sheet_path):
            raise FileNotFoundError(f"sampleSheet.tsv not found: {sample_sheet_path}")

        sample_order = []
        conditions = {}
        with open(sample_sheet_path, newline="") as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            for row in reader:
                sample = row["sampleName"].strip()
                cond   = row["condition"].strip()
                sample_order.append(sample)
                conditions[sample] = cond
        return sample_order, conditions

    # Main usage
    sample_sheet_path = os.path.join(input_dir, "sampleSheet.tsv")
    samples_detected = detect_samples(input_dir)
    sample_order, condition_map = load_conditions(sample_sheet_path)

    # Check for missing samples
    missing = [s for s in sample_order if s not in samples_detected]
    if missing:
        raise ValueError(f"Samples in sampleSheet but not in input_dir: {missing}")

    # Build ordered samples_dict with condition attached
    samples_dict = {
        s: {
            "R1": samples_detected[s]["R1"],
            "R2": samples_detected[s].get("R2"),
            "condition": condition_map[s]
        }
        for s in sample_order
    }


    # Print for verification
    print("Samples detected from directory & sampleSheet:")
    for sample, files in samples_dict.items():
        print(
            f"  {sample}: "
            f"R1 = {files['R1']}, "
            f"R2 = {files.get('R2', 'NA')}, "
            f"condition = {files['condition']}"
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
        f"pipeline={pipeline}",
        f"input_dir={input_dir}",
        f"output_dir={output_dir}",
        f"organism={organism}",
        f"genome_dir={genome_dir}",
        f"downsample_size={downsample_size}",
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
        "--cores", str(threads),
    ]

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
