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
    parser.add_argument("--restriction_enzyme", default="mboi", help="Restriction enzyme used in the HiChIP experiment. Options include 'mboi' (default), 'hindiii', 'dpnii', 'bglii', 'ncoi', 'msei', 'hinfI', 'mnase', 'arima'")
    parser.add_argument("--bin_size", type=int, default=5000, help="Bin size in base pairs for interaction analysis (default: 5000 bp).")
    parser.add_argument("--binning_range", type=int, default=2000000, help="Maximum distance in base pairs for binning interactions (default: 2000000 bp).")
    parser.add_argument("--length_cutoff", type=int, default=1000, help="Minimum fragment length in base pairs to include in analysis (default: 1000 bp).")
    parser.add_argument("--threads", type=int, default=30, help="Number of CPU threads to use for parallel processing (default: 30).")
    parser.add_argument("--fdr", type=float, default=0.01, help="Minimum FDR threshold for significant interactions (default: 0.01).")
    parser.add_argument("--macs2_peaks", default="NA", help="Optional path to a MACS2 peaks file (e.g., narrowPeak) from ChIP-seq data to refine interactions (default: none).")
    parser.add_argument("--fithichip_BiasType", default="1", help="FitHiChIP Bias correction type: coverage (1) or ICE (2) based. Default: 1.")
    parser.add_argument("--maps_count_cutoff", type=int, default=5, help="Minimum read count threshold for calling significant interactions in MAPS (default: 5).")
    parser.add_argument("--maps_ratio_cutoff", type=float, default=2.0, help="Minimum observed-to-expected ratio for significant interactions in MAPS (default: 2.0).")
    parser.add_argument("--maps_model", default="pospoisson", choices=["pospoisson", "negbinom"], help="Statistical model for regression analysis. Options: 'pospoisson' (default) or 'negbinom'.")
    parser.add_argument("--maps_sex_chroms", default="X", choices=["NA", "X", "Y", "XY"], help="Sex chromosomes to include in analysis. Options: 'NA' (none, default), 'X', 'Y', or 'XY'.")
    parser.add_argument("--hicpro_params", default="NA", help="Optional config file to pass to HiC-Pro. Default: NA.")
    parser.add_argument("--hichipper_params", default="NA", help="Optional space-separated parameters to pass to Hichipper (e.g., '--read-length 75'). Default: NA.")
    parser.add_argument("--hicdc_params", default="NA", help="Optional space-separated parameters to pass to Hicdcplus (e.g., '--PeakFile peaks.bed'). Default: NA.")
    parser.add_argument("--samples_comparison", default="NA NA", help="Sample names needed for differential analysis (e.g., --samples_comparison 'SampleA SampleB'). Default: NA.")
    parser.add_argument("--reference_condition", default="NA", help="Referece condition in differential analysis (e.g., --reference_condition 'Disease'). Default: NA.")
    parser.add_argument("--baseline_condition", default="NA", help="Baseline condition in differential analysis (e.g., --baseline_condition 'Healthy'). Default: NA.")
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

    

    restriction_enzyme=args.restriction_enzyme
    bin_size = args.bin_size
    binning_range = args.binning_range
    length_cutoff = args.length_cutoff
    threads = args.threads
    maps_model = args.maps_model
    maps_sex_chroms = args.maps_sex_chroms
    optical_duplicate_distance = 0  # Hardcoded default
    mapq = 30  # Hardcoded default
    generate_hic = 1  # Hardcoded default
    macs2_peaks = args.macs2_peaks
    hicpro_params = args.hicpro_params
    fithichip_BiasType = args.fithichip_BiasType
    hichipper_params = args.hichipper_params
    hicdc_params = args.hicdc_params
    samples_comparison = args.samples_comparison
    maps_count_cutoff = args.maps_count_cutoff
    maps_ratio_cutoff = args.maps_ratio_cutoff
    fdr = args.fdr
    method = args.method
    reference_condition = args.reference_condition
    baseline_condition = args.baseline_condition
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

    # Peak file detection
    if macs2_peaks != "NA":
        if not os.path.exists(macs2_peaks):
            print("Warning: The provided peaks file does not exist. Please provide a valid file path for called peaks.", file=sys.stderr)
            sys.exit(1)

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
        "--snakefile", os.path.join(script_dir, "../rules/Altpromoterflow.Snakefile"),
        "--printshellcmds",
        "--directory", output_dir,
        "--use-conda",
        "--conda-prefix", os.path.join(genome_dir, ".snakemake_conda"),
        #"--use-singularity",
        #"--singularity-args", f"-B {genome_dir} -B {singularity_tmpdir} -B {singularity_cachedir}",
        "--config",
        f"pipeline={pipeline}",
        f"input_dir={input_dir}",
        f"output_dir={output_dir}",
        f"organism={organism}",
        f"genome_dir={genome_dir}",
        f"downsample_size={downsample_size}",
        #f"samples={json.dumps(samples)}",
        f"reads={json.dumps(reads)}",
        f"restriction_enzyme={restriction_enzyme}",
        f"bin_size={bin_size}",
        f"binning_range={binning_range}",
        f"generate_hic={generate_hic}",
        f"mapq={mapq}",
        f"length_cutoff={length_cutoff}",
        f"threads={threads}",
        f"maps_model={maps_model}",
        f"maps_sex_chroms={maps_sex_chroms}",
        f"maps_count_cutoff={maps_count_cutoff}",
        f"maps_ratio_cutoff={maps_ratio_cutoff}",
        f"fdr={fdr}",
        f"optical_duplicate_distance={optical_duplicate_distance}",
        f"hicpro_params={hicpro_params}",
        f"fithichip_BiasType={fithichip_BiasType}",
        f"hichipper_params={hichipper_params}",
        f"hicdc_params={hicdc_params}",
        f"macs2_peaks={macs2_peaks}",
        f"samples_comparison={samples_comparison}",
        f"samples_dict={json.dumps(samples_dict)}",
        f"method={method}",
        #f"reads={','.join(reads)}",
        f"reference_condition={reference_condition}",
        f"baseline_condition={baseline_condition}",
        f"max_gFC={max_gFC}",
        f"min_pFC={min_pFC}",
        "--cores", str(threads),
    ]

    # Append any extra Snakemake arguments
    if extra_args:
        snakemake_cmd.extend(extra_args)

    # Generate unique log file name
    log_pattern = os.path.join(output_dir, "hichipsnake_run_*.log")
    log_number = max([int(re.search(r"hichipsnake_run_(\d+)\.log", f).group(1)) 
                    for f in glob.glob(log_pattern)], default=0) + 1
    log_file_path = os.path.join(output_dir, f"hichipsnake_run_{log_number}.log")

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
