import os
import re

print("CONFIG = ", config)
print("METHOD = ", config.get("method", "not provided"))

#############################################
# Configuration
#############################################
input_dir = config["input_dir"]
output_dir = config["output_dir"]
genome_dir = config["genome_dir"]
organism = config["organism"]
samples_dict = config["samples_dict"]  # dict: {sample: {R1: ..., R2: ..., condition: ...}}
samples = list(samples_dict.keys())
reads = config["reads"]  # ["R1", "R2"]
# BAM mode = all samples have a 'bam' key in samples_dict
bam_mode = all(("bam" in samples_dict[s]) for s in samples)
# 0 = off; 50/75/100/150 = simulate fixed read length
readlen_sim = int(config.get("readlen_sim", 0))

is_paired = "R2" in reads

if not bam_mode:
    if is_paired:
        ruleorder: dexseq_featurecounts_paired > dexseq_featurecounts_single
        ruleorder: trim_galore_paired > trim_galore_single
        ruleorder: star_paired > star_single
        ruleorder: salmon_quant_paired > salmon_quant_single
        ruleorder: simulate_readlen_paired > simulate_readlen_single
    else:
        ruleorder: dexseq_featurecounts_single > dexseq_featurecounts_paired
        ruleorder: trim_galore_single > trim_galore_paired
        ruleorder: star_single > star_paired
        ruleorder: salmon_quant_single > salmon_quant_paired


if is_paired:
    ruleorder: dexseq_featurecounts_paired > dexseq_featurecounts_single
    ruleorder: salmon_quant_paired > salmon_quant_single
else:
    ruleorder: dexseq_featurecounts_single > dexseq_featurecounts_paired
    ruleorder: salmon_quant_single > salmon_quant_paired

# General config
threads = config.get("threads", 16) / 7
downsample_size = config.get("downsample_size", 0)
fastqc = config.get("fastqc", False)
trim = config.get("trim", False)
trimmer_options = config.get("trimmer_options", "")
star_options = config.get("star_options", "")
test_condition = config.get("test_condition", "")
#batch = config.get("batch", "")
batch = config.get("all_batch_str", "")
control_condition = config.get("control_condition", "")
sheet_name = config.get("sheet_name", "")
max_gFC = config.get("max_gFC", 1.5)
min_pFC = config.get("min_pFC", 2.0)
lfcshrink = config.get("lfcshrink", False)
star_precomputed_dir = config.get("star_precomputed", "")

# Method selection
method = config.get("method", "rnaseq").lower()
do_salmon = method in ["salmon", "rnaseq"]
do_proactiv = method in ["proactiv", "rnaseq"]
do_dexseq = method in ["dexseq", "rnaseq"]
do_cage = method in ["cage"]  # only run cage if explicitly specified

# Group info
sample_conditions = [samples_dict[s]["condition"] for s in samples]

# Uniqued versions for comparisons
condition_compare = list(dict.fromkeys(sample_conditions))

# Construct strings
condition_str = ",".join(sample_conditions)
condition_compare_str = " ".join(condition_compare)

print("condition_str =", condition_str)
print("condition_compare_str =", condition_compare_str)

# batch = ["batch1"] * len(samples)
# batch_str = ",".join(batch)
batch_str = batch


def build_batch_condition(batch, samples=None):
    if batch is None:
        return ""
    # Convert batch to a list of cleaned tokens
    if isinstance(batch, str):
        s = batch.strip()
        if ",,," in s:
            parts = re.split(r",{3,}", s)  # split by 3 or more commas
            tokens = [p.replace(",", "").strip() for p in parts if p.strip()]
        else:
            tokens = [t.replace(",", "").strip() for t in re.split(r"[,\s]+", s) if t.strip()]
    else:
        # Already a list or tuple
        tokens = [str(t).strip() for t in batch if str(t).strip()]
    # Join all cleaned tokens into a single comma-separated string
    return ",".join(tokens)

batch_condition = build_batch_condition(batch, samples=samples)
print(batch_condition)
diff_samples_str    = config.get("diff_samples_str", "")
diff_samples = [s for s in config["diff_samples_str"].split(",") if s]
diff_conditions_str = config.get("diff_conditions_str", "")
diff_batch_str      = config.get("diff_batch_str", "")

diff_sample_conditions = [c for c in diff_conditions_str.split(",") if c]
unique_conditions = list(dict.fromkeys(diff_sample_conditions))
diff_condition_compare_str = " ".join(unique_conditions)

single_condition = (len(condition_compare) <= 1) or (test_condition == control_condition)
#do_diff = not single_condition  # run differential only when we truly have >=2 conditions
do_diff = bool(diff_samples_str.strip()) and (test_condition != control_condition)
print(f"do_diff (diff-based) = {do_diff}")
print(f"single_condition={single_condition}; do_diff={do_diff}")

comparison_name = f"{control_condition}_vs_{test_condition}"

def safe_name(s):
    return s.replace(" ", "_").replace("/", "_").replace(":", "_")

sheet_name = safe_name(sheet_name)
comparison_name = safe_name(comparison_name)

analysis_name = f"{sheet_name}__{comparison_name}"


# Paths from genomesetup.smk
star_index = f"{genome_dir}/organisms/{organism}/STARIndex"
genes_gtf = f"{genome_dir}/organisms/{organism}/Annotation.gtf"
genes_bed = f"{genome_dir}/organisms/{organism}/Annotation/genes.bed"
dexseq_gff = f"{genome_dir}/organisms/{organism}/Annotation/DEXSeq_flattened_exons.gff"
proactiv_rds = f"{genome_dir}/organisms/{organism}/Annotation/proActiv_promoter_annotation.rds"
tx2gene = f"{genome_dir}/organisms/{organism}/Annotation/genes_t2g.tsv"

# Path for sanity sanity_check
SCRIPTS_DIR = os.path.realpath(os.path.join(workflow.basedir, "../scripts"))

# Validate
if not all([output_dir, organism, samples]):
    raise ValueError("Missing required config variables: output_dir, organism, samples")
if not os.path.exists(genes_gtf):
    raise ValueError(f"GTF file {genes_gtf} does not exist; run genomesetup.smk first")

# Ensure output directories
os.makedirs(output_dir + "/fastqs/raw", exist_ok=True)
os.makedirs(output_dir + "/fastqs/downsampled", exist_ok=True)
os.makedirs(output_dir + "/fastqs/trimmed", exist_ok=True)
os.makedirs(output_dir + "/bam", exist_ok=True)
os.makedirs(output_dir + "/fastqc", exist_ok=True)
os.makedirs(output_dir + "/proactiv", exist_ok=True)
os.makedirs(output_dir + "/dexseq", exist_ok=True)
os.makedirs(output_dir + "/salmon", exist_ok=True)
os.makedirs("logs", exist_ok=True)
os.makedirs("benchmarks", exist_ok=True)

#############################################
# Final Target
#############################################
input_all = []

proactiv_merge_dir = f"{output_dir}/proactiv/counts_merged_diff/{analysis_name}"
dexseq_merge_dir= f"{output_dir}/dexseq/counts_merged_samplesheet/{analysis_name}"
salmon_merge_dir= f"{output_dir}/salmon/counts_merged_samplesheet/{analysis_name}"
cage_merge_dir = f"{output_dir}/cage/counts_merged_diff/{analysis_name}"

proactiv_merge_out = f"{output_dir}/proactiv/counts_merged_diff/{analysis_name}/proactiv_raw_counts.rds"
dexseq_merge_out = f"{output_dir}/dexseq/counts_merged_samplesheet/{analysis_name}/merged_promoter_counts.rds"
salmon_merge_out = f"{output_dir}/salmon/counts_merged_samplesheet/{analysis_name}/merged_promoter_counts.rds"
cage_merge_out = f"{output_dir}/cage/counts_merged_diff/{analysis_name}/cage_counts.rds"



proactiv_class_dir = f"{output_dir}/proactiv/promoter_classification_condition_wise/{analysis_name}"
dexseq_class_dir = f"{output_dir}/dexseq/promoter_classification_condition_wise/{analysis_name}"
salmon_class_dir = f"{output_dir}/salmon/promoter_classification_condition_wise/{analysis_name}"
cage_class_dir = f"{output_dir}/cage/promoter_classification_condition_wise/{analysis_name}"


def make_class_outputs(base_dir):
    return {
        "se":       f"{base_dir}/Promoter_activity_SE.rds",
        "rowData":  f"{base_dir}/Summary_classified_rowData.rds",
        "raw_count":f"{base_dir}/raw_promoter_counts.rds",
    }

proactiv_class_outputs = make_class_outputs(proactiv_class_dir)
dexseq_class_outputs   = make_class_outputs(dexseq_class_dir)
salmon_class_outputs   = make_class_outputs(salmon_class_dir)
cage_class_outputs   = make_class_outputs(cage_class_dir)



proactiv_diff_dir = f"{output_dir}/proactiv/differential/{analysis_name}"
dexseq_diff_dir = f"{output_dir}/dexseq/differential/{analysis_name}"
salmon_diff_dir = f"{output_dir}/salmon/differential/{analysis_name}"
cage_diff_dir = f"{output_dir}/cage/differential/{analysis_name}"

def make_diff_outputs(base_dir):
    return {
        "DEG": f"{base_dir}/Differential_promoter_usage_pFC2_gFC1_5.rds",
        "AP":  f"{base_dir}/Promoter_differential_activity_FDR0_05.rds",
    }

proactiv_diff_out = make_diff_outputs(proactiv_diff_dir)
dexseq_diff_out = make_diff_outputs(dexseq_diff_dir)
salmon_diff_out = make_diff_outputs(salmon_diff_dir)
cage_diff_out = make_diff_outputs(cage_diff_dir)




proactiv_plots_dir = f"{output_dir}/proactiv/plots_condition_wise/{analysis_name}"
dexseq_plots_dir   = f"{output_dir}/dexseq/plots_condition_wise/{analysis_name}"
salmon_plots_dir  = f"{output_dir}/salmon/plots_condition_wise/{analysis_name}"
cage_plots_dir  = f"{output_dir}/cage/plots_condition_wise/{analysis_name}"

def make_conditionwise_plot_outputs(plots_dir, unique_conditions):
    shared = expand(
        plots_dir + "/promoter_activity_{plot}.pdf",
        plot=["number_hist_all", "number_hist_without1", "tsne_plot"],
    )
    per_cond = [
        f"{plots_dir}/{cond}/{cond}_promoter_activity_{plot}.pdf"
        for cond in unique_conditions
        for plot in [
            "geneexpression_correlation", "category_percentage_genewise", "position_category",
            "category_percentage", "single_multiple_category", "category_comparison"
        ]
    ]
    return shared + per_cond

proactiv_plot_outputs = make_conditionwise_plot_outputs(proactiv_plots_dir, unique_conditions)
dexseq_plot_outputs   = make_conditionwise_plot_outputs(dexseq_plots_dir, unique_conditions)
salmon_plot_outputs   = make_conditionwise_plot_outputs(salmon_plots_dir, unique_conditions)
cage_plot_outputs   = make_conditionwise_plot_outputs(cage_plots_dir, unique_conditions)


if do_salmon:
    input_all += expand(output_dir + "/salmon/counts_separate/{sample}/quant.sf", sample=samples)
    input_all += expand(output_dir + "/salmon/counts_separate/{sample}_promoter_counts.rds", sample=samples)
    input_all.append(output_dir + "/salmon/counts_merged/merged_promoter_counts.rds")
    input_all.append(salmon_class_outputs["se"])
    input_all.append(output_dir + "/salmon/promoter_classification_total/Promoter_activity_SE.rds")
    if do_diff:
        input_all += list(salmon_diff_out.values())
    input_all += salmon_plot_outputs
    input_all += expand(output_dir + "/salmon/plots_total/overall/overall_promoter_activity_{plot}.pdf", 
                        plot=[
                            "category_percentage", "category_comparison", "position_category",
                            "geneexpression_correlation", "category_percentage_genewise",
                            "single_multiple_category", "number_hist_all", "number_hist_without1"
                        ])


if do_proactiv:
    input_all.append(output_dir + "/proactiv/counts_merged/proactiv_raw_counts.rds")
    input_all.append(proactiv_class_outputs["se"])
    input_all.append(output_dir + "/proactiv/promoter_classification_total/Promoter_activity_SE.rds")
    if do_diff:
        input_all += list(proactiv_diff_out.values())
    input_all += proactiv_plot_outputs
    input_all += expand(output_dir + "/proactiv/plots_total/overall/overall_promoter_activity_{plot}.pdf", 
                        plot=[
                            "category_percentage", "category_comparison", "position_category",
                            "geneexpression_correlation", "category_percentage_genewise",
                            "single_multiple_category", "number_hist_all", "number_hist_without1"
                        ])


if do_dexseq:
    input_all += expand(output_dir + "/STAR/{sample}/{sample}.strand.txt", sample=samples)
    input_all += expand(output_dir + "/dexseq/counts_separate/{sample}/{sample}_counts.txt", sample=samples)
    input_all += expand(output_dir + "/dexseq/counts_separate/{sample}/{sample}_promoter_counts.rds", sample=samples)
    input_all.append(output_dir + "/dexseq/counts_merged/merged_promoter_counts.rds")
    input_all.append(dexseq_class_outputs["se"])
    input_all.append(output_dir + "/dexseq/promoter_classification_total/Promoter_activity_SE.rds")
    if do_diff:
        input_all += list(dexseq_diff_out.values())
    input_all += dexseq_plot_outputs
    input_all += expand(output_dir + "/dexseq/plots_total/overall/overall_promoter_activity_{plot}.pdf", 
                        plot=[
                            "category_percentage", "category_comparison", "position_category",
                            "geneexpression_correlation", "category_percentage_genewise",
                            "single_multiple_category", "number_hist_all", "number_hist_without1"
                        ])



if do_cage:
    input_all += expand(output_dir + "/cage/counts_separate/{sample}_promoter_counts.txt", sample=samples)
    input_all.append(output_dir + "/cage/counts_merged/cage_counts.rds")
    input_all.append(cage_class_outputs["se"])
    input_all.append(output_dir + "/cage/promoter_classification_total/Promoter_activity_SE.rds")
    if do_diff:
        input_all += list(cage_diff_out.values())
    input_all += cage_plot_outputs
    input_all += expand(output_dir + "/cage/plots_total/overall/overall_promoter_activity_{plot}.pdf", 
                        plot=[
                            "category_percentage", "category_comparison", "position_category",
                            "geneexpression_correlation", "category_percentage_genewise",
                            "single_multiple_category", "number_hist_all", "number_hist_without1"
                        ])


input_all.append(output_dir + "/multiqc/multiqc_report.html")

rule all:   
    input:
        #expand(output_dir + "/fastqs/raw/{sample}/{sample}_{read}.fastq.gz", sample=samples, read=reads),
        expand(output_dir + "/STAR/{sample}/{sample}.sorted.bam", sample=samples),
        expand(output_dir + "/bigwig/{sample}.bw", sample=samples),
        input_all


#############################################
# Preprocessing Rules
#############################################
if not bam_mode:
    rule link_fastqs:
        input:
            fastqs = lambda wildcards: config["samples_dict"][wildcards.sample][wildcards.read]
        output:
            fastq = output_dir + "/fastqs/raw/{sample}/{sample}_{read}.fastq.gz"
        log:
            output_dir + "/logs/link_fastqs.{sample}_{read}.log"
        benchmark:
            "benchmarks/link_fastqs.{sample}_{read}.txt"
        conda: "envs/altbasic.yaml"
        shell:
            """
            mkdir -p $(dirname {output.fastq})
            if [ ! -s "{input.fastqs}" ]; then
                echo "[WARN] No input FASTQ for {wildcards.sample} {wildcards.read} ({input.fastqs}); creating empty placeholder and skipping link." > {log}
                # empty file
                touch {output.fastq}
                exit 0
            fi
            ln -sf "{input.fastqs}" "{output.fastq}" 2> "{log}"
            """

    rule fastqc:
        input:
            fastq=lambda wildcards: (
                output_dir + "/fastqs/downsampled/" + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.read + ".fastq.gz"
                if downsample_size > 0
                else output_dir + "/fastqs/raw/" + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.read + ".fastq.gz"
            )
        output:
            html=output_dir + "/fastqc/{sample}/{sample}_{read}_fastqc.html",
            zip=output_dir + "/fastqc/{sample}/{sample}_{read}_fastqc.zip"
        log:
            "logs/fastqc.{sample}_{read}.log"
        threads: threads
        benchmark:
            "benchmarks/fastqc.{sample}_{read}.txt"
        conda: "envs/altbasic.yaml"
        shell:
            """
            mkdir -p {output_dir}/fastqc/{wildcards.sample}
            if [ ! -s {input.fastq} ]; then
                echo "[WARN] Empty/missing input FASTQ {input.fastq}; skipping FastQC for {wildcards.sample}_{wildcards.read}." > {log}
                touch {output.html} {output.zip}
                exit 0
            fi

            base=$(basename {input.fastq} .fastq.gz)
            fastqc -t {threads} -o {output_dir}/fastqc {input.fastq} > {log} 2>&1
            mv {output_dir}/fastqc/${{base}}_fastqc.html {output.html}
            mv {output_dir}/fastqc/${{base}}_fastqc.zip {output.zip}
            """

    rule downsample:
        input:
            fastq=output_dir + "/fastqs/raw/{sample}/{sample}_{read}.fastq.gz"
        output:
            fastq=output_dir + "/fastqs/downsampled/{sample}/{sample}_{read}.fastq.gz"
        params:
            size=downsample_size
        log:
            "logs/downsample.{sample}_{read}.log"
        benchmark:
            "benchmarks/downsample.{sample}_{read}.txt"
        threads: threads
        conda: "envs/altbasic.yaml"
        shell:
            """
            mkdir -p $(dirname {output.fastq})
            if [ ! -s {input.fastq} ]; then
                echo "[WARN] Empty/missing input FASTQ {input.fastq}; skipping downsample and creating empty placeholder." > {log}
                touch {output.fastq}
                exit 0
            fi
            if [ {params.size} -gt 0 ]; then
                seqtk sample -s100 {input.fastq} {params.size} | gzip > {output.fastq} 2> {log}
            else
                ln -sf {input.fastq} {output.fastq} 2> {log}
            fi
            """

    rule trim_galore_single:
        input:
            r1 = output_dir + "/fastqs/downsampled/{sample}/{sample}_R1.fastq.gz"
        output:
            r1        = output_dir + "/fastqs/trimmed/{sample}/{sample}_R1.fastq.gz",
            report_r1 = output_dir + "/fastqs/trimmed/{sample}/{sample}_R1_fastqc.html"
        params:
            opts   = trimmer_options,
            outdir = output_dir + "/fastqs/trimmed"
        log:
            "logs/trim_galore.{sample}.log"
        benchmark:
            "benchmarks/trim_galore.{sample}.txt"
        threads: threads
        conda: "envs/altbasic.yaml"
        shell:
            """
            mkdir -p $(dirname {output.r1})

            if [ ! -s {input.r1} ]; then
                echo "[WARN] Empty/missing input FASTQ {input.r1}; skipping trim_galore (single-end) for {wildcards.sample}." > {log}
                # Placeholder
                touch {output.r1}
                touch {output.report_r1}
                exit 0
            fi

            trim_galore --stringency 3 --cores {threads} {params.opts} \
                --output_dir {params.outdir} --gzip \
                --fastqc --fastqc_args "-o {params.outdir}" \
                {input.r1} > {log} 2>&1

            mv {params.outdir}/{wildcards.sample}_R1_trimmed.fq.gz  {output.r1}
            mv {params.outdir}/{wildcards.sample}_R1_trimmed_fastqc.html     {output.report_r1}
            """

    rule trim_galore_paired:
        input:
            r1=output_dir + "/fastqs/downsampled/{sample}/{sample}_R1.fastq.gz",
            r2=output_dir + "/fastqs/downsampled/{sample}/{sample}_R2.fastq.gz"
        output:
            r1=output_dir + "/fastqs/trimmed/{sample}/{sample}_R1.fastq.gz",
            r2=output_dir + "/fastqs/trimmed/{sample}/{sample}_R2.fastq.gz",
            report_r1=output_dir + "/fastqs/trimmed/{sample}/{sample}_R1_fastqc.html",
            report_r2=output_dir + "/fastqs/trimmed/{sample}/{sample}_R2_fastqc.html"
        params:
            opts=trimmer_options,
            outdir=output_dir + "/fastqs/trimmed"
        log:
            "logs/trim_galore.{sample}.log"
        benchmark:
            "benchmarks/trim_galore.{sample}.txt"
        threads: threads
        conda: "envs/altbasic.yaml"
        shell:
            """
            mkdir -p $(dirname {output.r1}) $(dirname {output.r2})

            if [ ! -s {input.r1} ] || [ ! -s {input.r2} ]; then
                echo "[WARN] Empty/missing input FASTQs {input.r1} / {input.r2}; skipping trim_galore (paired-end) for {wildcards.sample}." > {log}
                touch {output.r1} {output.r2} {output.report_r1} {output.report_r2}
                exit 0
            fi

            trim_galore --paired --stringency 3 --cores {threads} {params.opts} \
                --output_dir {params.outdir} --gzip \
                --fastqc --fastqc_args "-o {params.outdir}" \
                {input.r1} {input.r2} > {log} 2>&1
            mv {params.outdir}/{wildcards.sample}_R1_val_1.fq.gz {output.r1}
            mv {params.outdir}/{wildcards.sample}_R2_val_2.fq.gz {output.r2}
            mv {params.outdir}/{wildcards.sample}_R1_val_1_fastqc.html {output.report_r1}
            mv {params.outdir}/{wildcards.sample}_R2_val_2_fastqc.html {output.report_r2}
            """

if bam_mode:
    rule star_import_bam:
        input:
            bam = lambda wc: samples_dict[wc.sample]["bam"],
            sj  = lambda wc: samples_dict[wc.sample]["sj"]
        output:
            bam = output_dir + "/STAR/{sample}/{sample}.sorted.bam",
            sj  = output_dir + "/STAR/junctions/{sample}.SJ.out.tab"
        log:
            "logs/star_import_bam.{sample}.log"
        conda:
            "envs/bam.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output.bam}) $(dirname {output.sj})
            echo "[INFO] Importing precomputed BAM and SJ for {wildcards.sample}" > {log}
            ln -sf {input.bam} {output.bam}
            ln -sf {input.sj}  {output.sj}
            """

    rule star_import_bai:
            input:
                bai = lambda wc: samples_dict[wc.sample]["bai"]
            output:
                bai = output_dir + "/STAR/{sample}/{sample}.sorted.bam.bai",
            log:
                "logs/star_import_bai.{sample}.log"
            conda:
                "envs/bam.yaml"
            shell:
                r"""
                mkdir -p $(dirname {output.bai})
                echo "[INFO] Importing precomputed BAI for {wildcards.sample}" > {log}
                ln -sf {input.bai} {output.bai}
                """

else:
    rule simulate_readlen_single:
        input:
            r1 = lambda wc: (
                output_dir + f"/fastqs/trimmed/{wc.sample}/{wc.sample}_R1.fastq.gz"
                if trim else
                output_dir + f"/fastqs/downsampled/{wc.sample}/{wc.sample}_R1.fastq.gz"
            )
        output:
            r1 = output_dir + f"/fastqs/readlen_{readlen_sim}" + "/{sample}/{sample}_R1.fastq.gz"
        params:
            L = readlen_sim
        log:
            "logs/readlen_sim.{sample}.log"
        threads: 1
        conda: "envs/altbasic.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.r1})
            if [ {params.L} -le 0 ]; then
                ln -sf {input.r1} {output.r1}
                exit 0
            fi
            seqtk trimfq -L {params.L} {input.r1} | gzip > {output.r1}
            """

    rule star_single:
        input:
            r1=lambda wc: (
                output_dir + f"/fastqs/readlen_{readlen_sim}/{wc.sample}/{wc.sample}_R1.fastq.gz"
                if readlen_sim > 0 else
                (
                    output_dir + f"/fastqs/trimmed/{wc.sample}/{wc.sample}_R1.fastq.gz"
                    if trim else
                    output_dir + f"/fastqs/downsampled/{wc.sample}/{wc.sample}_R1.fastq.gz"
                )
            ),
            index=star_index,
            gtf=genes_gtf
        output:
            bam=temp(output_dir + "/STAR/{sample}/{sample}.sorted.bam"),
            sj=output_dir + "/STAR/junctions/{sample}.SJ.out.tab"
        params:
            opts=star_options,
            prefix=output_dir + "/STAR/junctions/{sample}.",
            sample_dir=output_dir + "/STAR/{sample}",
            samsort_memory="8G",
            samtools_threads=16,
            precomputed_dir=star_precomputed_dir
        log:
            "logs/star.{sample}.log"
        benchmark:
            "benchmarks/star.{sample}.txt"
        threads: min(32, threads)
        conda: "envs/basic.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p {params.sample_dir}
            mkdir -p $(dirname {output.sj})

            STAR --runThreadN {threads} \
                {params.opts} \
                --sjdbOverhang 100 \
                --outSAMunmapped Within \
                --outSAMtype BAM Unsorted \
                --sjdbGTFfile {input.gtf} \
                --genomeDir {input.index} \
                --readFilesIn {input.r1} \
                --readFilesCommand 'gunzip -c' \
                --outFileNamePrefix {params.prefix} >> {log} 2>&1

            samtools sort -m {params.samsort_memory} \
                -T {params.sample_dir}/{wildcards.sample} \
                -@ {params.samtools_threads} \
                -O bam \
                -o {output.bam} \
                {params.prefix}Aligned.out.bam >> {log} 2>&1
            """


    rule simulate_readlen_paired:
        input:
            r1=lambda wc: (
                output_dir + f"/fastqs/trimmed/{wc.sample}/{wc.sample}_R1.fastq.gz"
                if trim else
                output_dir + f"/fastqs/downsampled/{wc.sample}/{wc.sample}_R1.fastq.gz"
            ),
            r2=lambda wc: (
                output_dir + f"/fastqs/trimmed/{wc.sample}/{wc.sample}_R2.fastq.gz"
                if trim else
                output_dir + f"/fastqs/downsampled/{wc.sample}/{wc.sample}_R2.fastq.gz"
            )
        output:
            r1=output_dir + f"/fastqs/readlen_{readlen_sim}" + "/{sample}/{sample}_R1.fastq.gz",
            r2=output_dir + f"/fastqs/readlen_{readlen_sim}" + "/{sample}/{sample}_R2.fastq.gz"
        params:
            L=readlen_sim
        log:
            "logs/readlen_sim.{sample}.log"
        threads: 1
        conda: "envs/altbasic.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.r1})

            if [ {params.L} -le 0 ]; then
                ln -sf {input.r1} {output.r1}
                ln -sf {input.r2} {output.r2}
                exit 0
            fi

            seqtk trimfq -L {params.L} {input.r1} | gzip > {output.r1}
            seqtk trimfq -L {params.L} {input.r2} | gzip > {output.r2}
            """


    rule star_paired:
        input:
            r1=lambda wc: (
                output_dir + f"/fastqs/readlen_{readlen_sim}/{wc.sample}/{wc.sample}_R1.fastq.gz"
                if readlen_sim > 0 else
                (
                    output_dir + f"/fastqs/trimmed/{wc.sample}/{wc.sample}_R1.fastq.gz"
                    if trim else
                    output_dir + f"/fastqs/downsampled/{wc.sample}/{wc.sample}_R1.fastq.gz"
                )
            ),
            r2=lambda wc: (
                output_dir + f"/fastqs/readlen_{readlen_sim}/{wc.sample}/{wc.sample}_R2.fastq.gz"
                if readlen_sim > 0 else
                (
                    output_dir + f"/fastqs/trimmed/{wc.sample}/{wc.sample}_R2.fastq.gz"
                    if trim else
                    output_dir + f"/fastqs/downsampled/{wc.sample}/{wc.sample}_R2.fastq.gz"
                )
            ),
            index=star_index,
            gtf=genes_gtf
        output:
            bam=temp(output_dir + "/STAR/{sample}/{sample}.sorted.bam"),
            sj=output_dir + "/STAR/junctions/{sample}.SJ.out.tab"
        params:
            opts=star_options,
            prefix=output_dir + "/STAR/junctions/{sample}.",
            sample_dir=output_dir + "/STAR/{sample}",
            samsort_memory="8G",
            samtools_threads=16,
            precomputed_dir=star_precomputed_dir
        log:
            "logs/star.{sample}.log"
        benchmark:
            "benchmarks/star.{sample}.txt"
        threads: min(32, threads)
        conda: "envs/basic.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p {params.sample_dir}
            mkdir -p $(dirname {output.sj})

            STAR --runThreadN {threads} \
                {params.opts} \
                --sjdbOverhang 100 \
                --outSAMunmapped Within \
                --outSAMtype BAM Unsorted \
                --sjdbGTFfile {input.gtf} \
                --genomeDir {input.index} \
                --readFilesIn {input.r1} {input.r2} \
                --readFilesCommand 'gunzip -c' \
                --outFileNamePrefix {params.prefix} >> {log} 2>&1

            samtools sort -m {params.samsort_memory} \
                -T {params.sample_dir}/{wildcards.sample} \
                -@ {params.samtools_threads} \
                -O bam \
                -o {output.bam} \
                {params.prefix}Aligned.out.bam >> {log} 2>&1
            """
    rule index_bam:
        input:
            bam = output_dir + "/STAR/{sample}/{sample}.sorted.bam"
        output:
            bai = output_dir + "/STAR/{sample}/{sample}.sorted.bam.bai"
        params:
            precomputed_dir = star_precomputed_dir
        conda: "envs/bam.yaml"
        log:
            "logs/index_bam.{sample}.log"
        benchmark:
            "benchmarks/index_bam.{sample}.txt"
        threads: 4
        shell:
            """
            if [ -n "{params.precomputed_dir}" ] && \
                    [ -s "{params.precomputed_dir}/{wildcards.sample}.sorted.bam.bai" ]; then
                        echo "Using precomputed BAI for {wildcards.sample}" >> {log}
                        ln -sf {params.precomputed_dir}/{wildcards.sample}.sorted.bam.bai {output.bai}
                        exit 0
            fi
            samtools index -@ {threads} {input.bam} 2> {log}
            """

rule bam_to_bigwig:
    input:
        bam = output_dir + "/STAR/{sample}/{sample}.sorted.bam",
        bai = output_dir + "/STAR/{sample}/{sample}.sorted.bam.bai"
    output:
        bigwig = output_dir + "/bigwig/{sample}.bw"
    params:
        binsize = 50
    conda: "envs/bam.yaml"
    log:
        "logs/bam_to_bigwig.{sample}.log"
    benchmark:
        "benchmarks/bam_to_bigwig.{sample}.txt"
    threads: threads
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} \
            --binSize {params.binsize} --normalizeUsing RPKM \
            --numberOfProcessors {threads} 2> {log}
        """

#############################################
# proActiv Rules
#############################################

rule proactiv_count:
    output:
        junctions = output_dir + "/proactiv/counts_merged/proactiv_raw_counts.rds"
    input:
        promoter_rds = proactiv_rds,
        sj_files = [f"{output_dir}/STAR/junctions/{s}.SJ.out.tab" for s in samples]
    params:
        condition = condition_str,
        sj_list   = output_dir + "/proactiv/counts_merged/sj_files.list"
    log:
        "logs/proactiv_counts.log"
    threads: threads
    conda: "envs/dexR.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {output_dir}/proactiv/counts_merged

        echo "[INFO] writing SJ list" | tee {log}

        printf "%s\n" {input.sj_files} > {params.sj_list}

        echo "[INFO] SJ list written to {params.sj_list}" | tee -a {log}
        head -n 3 {params.sj_list} | tee -a {log}

        echo "[INFO] running Rscript..." | tee -a {log}
        Rscript {workflow.basedir}/../scripts/proactiv_count_merged.R \
            {output_dir}/proactiv/counts_merged \
            {input.promoter_rds} \
            {params.sj_list} \
            "{params.condition}" 2>&1 | tee -a {log}

        echo "[INFO] done" | tee -a {log}
        """


rule proactiv_count_diff:
    output:
        junctions = proactiv_merge_out
    input:
        promoter_rds = proactiv_rds,
        sj_files = [f"{output_dir}/STAR/junctions/{s}.SJ.out.tab" for s in diff_samples]
    params:
        out_dir = proactiv_merge_dir,
        condition = diff_conditions_str,
        sj_list   = output_dir + "/proactiv/counts_merged_diff/sj_files.list"
    log:
        "logs/proactiv_counts.log"
    threads: threads
    conda: "envs/dexR.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.out_dir}"

        echo "[INFO] writing SJ list" | tee {log}

        printf "%s\n" {input.sj_files} > {params.sj_list}

        echo "[INFO] SJ list written to {params.sj_list}" | tee -a {log}
        head -n 3 {params.sj_list} | tee -a {log}

        echo "[INFO] running Rscript..." | tee -a {log}
        Rscript {workflow.basedir}/../scripts/proactiv_count_merged.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            {params.sj_list} \
            "{params.condition}" 2>&1 | tee -a {log}

        echo "[INFO] done" | tee -a {log}
        """


rule proactiv_promoter_classification_total:
    output:
        rowData = output_dir + "/proactiv/promoter_classification_total/Summary_classified_rowData.rds",
        SE=output_dir + "/proactiv/promoter_classification_total/Promoter_activity_SE.rds"
    input:
        combined=output_dir + "/proactiv/counts_merged/proactiv_raw_counts.rds",
        promoter_rds=proactiv_rds
    params:
        samples=" ".join(samples),
        newnames=" ".join(samples),
        condition = condition_str,
        condition_compare = condition_compare_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        batch=batch_condition,
        norm_method=config.get("norm_method", "deseq2")

    log:
        "logs/proactiv_promoter_classification_total.log"
    benchmark:
        "benchmarks/promoter_classification_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/promoter_classify_total.R \
            {output_dir}/proactiv/promoter_classification_total \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.combined} \
            "{params.condition_compare}" \
            "{params.fit_script}" \
            "{params.batch}" \
            "{params.norm_method}" > {log} 2>&1
        """

rule proactiv_promoter_classification:
    output:
        se       = proactiv_class_outputs["se"],
        rowData  = proactiv_class_outputs["rowData"],
        raw_count= proactiv_class_outputs["raw_count"]
    input:
        promoter_rds=proactiv_rds,
        junctions=proactiv_merge_out,
    params:
        out_dir          = proactiv_class_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        batch = diff_batch_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/proactiv_promoter_classification_condition_wise.log"
    benchmark:
        "benchmarks/proactiv_promoter_classification_condition_wise.txt"
    threads: threads
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_classify.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.junctions} \
            "{params.condition_compare}" \
            "{params.fit_script}" \
            "{params.batch}" \
            "{params.norm_method}" > {log} 2>&1
        """

rule proactiv_differential:
    input:
        promoter_rds  = proactiv_rds,
        raw_count= proactiv_class_outputs["raw_count"],
        SE=proactiv_class_outputs["se"]
    output:
        deg = proactiv_diff_out["DEG"],
        altpromoter    = proactiv_diff_out["AP"]
    params:
        out_dir = proactiv_diff_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        baseline      = control_condition,
        reference    = test_condition,
        min_promoter_fold_change = min_pFC,
        max_gene_fold_change = max_gFC,
        lfcshrink = lfcshrink,
        batch = diff_batch_str,
    log:
        "logs/proactiv_promoter_differential.log"
    benchmark:
        "benchmarks/proactiv_promoter_differential.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_differential.R \
           "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.raw_count} \
            {input.SE} \
            {params.baseline} \
            {params.reference} \
            {params.min_promoter_fold_change} \
            {params.max_gene_fold_change} {params.lfcshrink} {params.batch} > {log} 2>&1 
        """

rule proactiv_promoter_plots:
    input:
        rowData      = proactiv_class_outputs["rowData"],
        promoter_rds = proactiv_rds,
        SE           = proactiv_class_outputs["se"]
    output:
        plots = proactiv_plot_outputs
    params:
        out_dir = proactiv_plots_dir,
        condition_compare = diff_condition_compare_str,
        condition = diff_conditions_str,
    log:
        "logs/proactiv_plots.log"
    benchmark:
        "benchmarks/proactiv_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/plots.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """

proactiv_overall_plot_outputs = [
    f"{output_dir}/proactiv/plots_total/overall/overall_promoter_activity_{plot}.pdf"
    for plot in [
        "category_percentage", "category_comparison", "position_category",
        "geneexpression_correlation", "category_percentage_genewise",
        "single_multiple_category",  # Plot 1–6
        "number_hist_all", "number_hist_without1"  # Plot 7–9
    ]
]

rule proactiv_overall_plots:
    input:
        rowData = output_dir + "/proactiv/promoter_classification_total/Summary_classified_rowData.rds",
        promoter_rds = proactiv_rds,
        SE=output_dir + "/proactiv/promoter_classification_total/Promoter_activity_SE.rds"
    output:
        plots = proactiv_overall_plot_outputs
    params:
        condition = condition_str,
        condition_compare = condition_compare_str,
    log:
        "logs/proactiv_plots_total.log"
    benchmark:
        "benchmarks/proactiv_plots_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p {output_dir}/proactiv/plots_total

        Rscript {workflow.basedir}/../scripts/plots_total.R \
            {output_dir}/proactiv/plots_total \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """

#############################################
# DEXSeq Rules
#############################################
rule infer_strand:
    input:
        bam=output_dir + "/STAR/{sample}/{sample}.sorted.bam",
        bed=genes_bed
    output:
        strand=output_dir + "/STAR/{sample}/{sample}.strand.txt"
    log:
        "logs/infer_strand.{sample}.log"
    benchmark:
        "benchmarks/infer_strand.{sample}.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        infer_experiment.py -r {input.bed} -i {input.bam} > {output.strand} 2> {log}
        """

rule dexseq_featurecounts_single:
    input:
        bam=output_dir + "/STAR/{sample}/{sample}.sorted.bam",
        strand=output_dir + "/STAR/{sample}/{sample}.strand.txt",
        gff=dexseq_gff
    output:
        counts=output_dir + "/dexseq/counts_separate/{sample}/{sample}_counts.txt"
    params:
        threads=16
    log:
        "logs/dexseq_featurecounts.{sample}.log"
    benchmark:
        "benchmarks/dexseq_featurecounts.{sample}.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p $(dirname {output.counts})

        strand=$(python {workflow.basedir}/../scripts/get_strand.py {input.strand})
        featureCounts -f -O -s $strand -T {params.threads} -F GFF \
            -t exonic_part -g gene_id \
            -a {input.gff} -o {output.counts} {input.bam} > {log} 2>&1
        """

rule dexseq_featurecounts_paired:
    input:
        bam=output_dir + "/STAR/{sample}/{sample}.sorted.bam",
        strand=output_dir + "/STAR/{sample}/{sample}.strand.txt",
        gff=dexseq_gff
    output:
        counts=output_dir + "/dexseq/counts_separate/{sample}/{sample}_counts.txt"
    params:
        threads=16
    log:
        "logs/dexseq_featurecounts.{sample}.log"
    benchmark:
        "benchmarks/dexseq_featurecounts.{sample}.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p $(dirname {output.counts})

        strand=$(python {workflow.basedir}/../scripts/get_strand.py {input.strand})
        featureCounts -f -O -s $strand -p -T {params.threads} -F GFF \
            -t exonic_part -g gene_id \
            -a {input.gff} -o {output.counts} {input.bam} > {log} 2>&1
        """

rule dexseq_counts_separate:
    output:
        counts=output_dir + "/dexseq/counts_separate/{sample}/{sample}_promoter_counts.rds"
    input:
        gff=dexseq_gff,
        counts=output_dir + "/dexseq/counts_separate/{sample}/{sample}_counts.txt",
        promoter_rds=proactiv_rds
    params:
        samples="{sample}",
        newnames="{sample}"
    log:
        "logs/dexseq_counts_separate.{sample}.log"
    benchmark:
        "benchmarks/dexseq_counts_separate.{sample}.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/dexseq_counts_separate.R \
            {output_dir}/dexseq/counts_separate \
            {input.gff} \
            {input.counts} \
            {input.promoter_rds} \
            {params.samples} \
            {params.newnames} > {log} 2>&1
        """


rule dexseq_promoter_merge:
    input:
        counts=expand(output_dir + "/dexseq/counts_separate/{sample}/{sample}_promoter_counts.rds", sample=samples),
        promoter_rds = proactiv_rds
    output:
        combined=output_dir + "/dexseq/counts_merged/merged_promoter_counts.rds"
    params:
        samples=samples,
        condition = condition_str,
        condition_compare = condition_compare_str,
        newnames=samples,
    log:
          "logs/dexseq_counts_merge.log"    
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/merge_promoter_counts_copy.R \
            {output_dir}/dexseq/counts_merged \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            "{params.condition_compare}"\
            {input.counts}  > {log} 2>&1
        """

rule dexseq_promoter_merge_samplesheet:
    input:
        counts=expand(output_dir + "/dexseq/counts_separate/{sample}/{sample}_promoter_counts.rds", sample = diff_samples),
        promoter_rds = proactiv_rds
    output:
        combined = dexseq_merge_out
    params:
        out_dir = dexseq_merge_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        newnames = diff_samples_str,
    log:
          "logs/dexseq_counts_merge.log"    
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"
        Rscript {workflow.basedir}/../scripts/merge_promoter_counts.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            "{params.condition_compare}"\
            {input.counts}  > {log} 2>&1
        """

rule dexseq_promoter_classification_total:
    output:
        rowData = output_dir + "/dexseq/promoter_classification_total/Summary_classified_rowData.rds",
        SE=output_dir + "/dexseq/promoter_classification_total/Promoter_activity_SE.rds"
    input:
        combined=output_dir + "/dexseq/counts_merged/merged_promoter_counts.rds",
        promoter_rds=proactiv_rds
    params:
        samples=" ".join(samples),
        newnames=" ".join(samples),
        condition = condition_str,
        condition_compare = condition_compare_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        batch=batch_condition,
        norm_method=config.get("norm_method", "deseq2")

    log:
        "logs/dexseq_promoter_classification_total.log"
    benchmark:
        "benchmarks/dexseq_promoter_classification_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/promoter_classify_total.R \
            {output_dir}/dexseq/promoter_classification_total \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.combined} \
            "{params.condition_compare}" \
            "{params.fit_script}" \
            "{params.batch}" \
            "{params.norm_method}" > {log} 2>&1
        """

rule dexseq_promoter_classification:
    input:
        combined = dexseq_merge_out,
        promoter_rds=proactiv_rds
    output:
        se       = dexseq_class_outputs["se"],
        rowData  = dexseq_class_outputs["rowData"],
        raw_count= dexseq_class_outputs["raw_count"]
    params:
        out_dir = dexseq_class_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        batch_unsorted = diff_batch_str,
        newnames=samples,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        test_condition = test_condition,
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/dexseq_promoter_classification.log"
    benchmark:
        "benchmarks/dexseq_promoter_classification.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_classify.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.combined} \
            "{params.condition_compare}"\
            "{params.fit_script}" \
            "{params.batch_unsorted}" \
            {params.norm_method} > {log} 2>&1
        """

rule dexseq_differential:
    input:
        promoter_rds  = proactiv_rds,
        raw_count= dexseq_class_outputs["raw_count"],
        SE=dexseq_class_outputs["se"]
    output:
        deg = dexseq_diff_out["DEG"],
        altpromoter    = dexseq_diff_out["AP"]
    params:
        out_dir = dexseq_diff_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        baseline      = control_condition,
        reference    = test_condition,
        min_promoter_fold_change = min_pFC,
        max_gene_fold_change = max_gFC,
        lfcshrink = lfcshrink,
        batch = diff_batch_str,
    log:
        "logs/dexseq_promoter_differential.log"
    benchmark:
        "benchmarks/dexseq_promoter_differential.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_differential.R \
           "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.raw_count} \
            {input.SE} \
            {params.baseline} \
            {params.reference} \
            {params.min_promoter_fold_change} \
            {params.max_gene_fold_change} {params.lfcshrink} {params.batch} > {log} 2>&1 
        """

rule dexseq_plots:
    input:
        rowData      = dexseq_class_outputs["rowData"],
        promoter_rds = proactiv_rds,
        SE           = dexseq_class_outputs["se"]
    output:
        plots = dexseq_plot_outputs
    params:
        out_dir = dexseq_plots_dir,
        condition_compare = diff_condition_compare_str,
        condition = diff_conditions_str,
    log:
        "logs/dexseq_plots.log"
    benchmark:
        "benchmarks/dexseq_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/plots.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """

dexseq_overall_plot_outputs = [
    f"{output_dir}/dexseq/plots_total/overall/overall_promoter_activity_{plot}.pdf"
    for plot in [
        "category_percentage", "category_comparison", "position_category",
        "geneexpression_correlation", "category_percentage_genewise",
        "single_multiple_category",  # Plot 1–6
        "number_hist_all", "number_hist_without1"  # Plot 7–9
    ]
]

rule dexseq_overall_plots:
    input:
        rowData = output_dir + "/dexseq/promoter_classification_total/Summary_classified_rowData.rds",
        promoter_rds = proactiv_rds,
        SE=output_dir + "/dexseq/promoter_classification_total/Promoter_activity_SE.rds"
    output:
        plots = dexseq_overall_plot_outputs
    params:
        condition = condition_str,
        condition_compare = condition_compare_str,
    log:
        "logs/dexseq_plots_total.log"
    benchmark:
        "benchmarks/dexseq_plots_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p {output_dir}/dexseq/plots_total

        Rscript {workflow.basedir}/../scripts/plots_total.R \
            {output_dir}/dexseq/plots_total \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """

#############################################
# Salmon Rules
#############################################
rule salmon_quant_single:
    input:
        r1=lambda wildcards: (
            output_dir + f"/fastqs/trimmed/{wildcards.sample}/{wildcards.sample}_R1.fastq.gz"
            if trim
            else output_dir + f"/fastqs/downsampled/{wildcards.sample}/{wildcards.sample}_R1.fastq.gz"
        )
    output:
        quant=output_dir + "/salmon/counts_separate/{sample}/quant.sf"
    params:
        outdir=output_dir + "/salmon/counts_separate/{sample}",
        libtype="A",  # Automatic library type detection
        index=f"{genome_dir}/organisms/{organism}/SalmonIndex"
    log:
        "logs/salmon_quant.{sample}.log"
    benchmark:
        "benchmarks/salmon_quant.{sample}.txt"
    threads: threads
    conda: "envs/basic.yaml"
    shell:
        """
        salmon quant -i {params.index} -l {params.libtype} \
            -r {input.r1} \
            -p {threads} -o {params.outdir} --gcBias --validateMappings > {log} 2>&1
        """

rule salmon_quant_paired:
    input:
        r1=lambda wildcards: (
            output_dir + f"/fastqs/trimmed/{wildcards.sample}/{wildcards.sample}_R1.fastq.gz"
            if trim
            else output_dir + f"/fastqs/downsampled/{wildcards.sample}/{wildcards.sample}_R1.fastq.gz"
        ),
        r2=lambda wildcards: (
            output_dir + f"/fastqs/trimmed/{wildcards.sample}/{wildcards.sample}_R2.fastq.gz"
            if trim
            else output_dir + f"/fastqs/downsampled/{wildcards.sample}/{wildcards.sample}_R2.fastq.gz"
        )
    output:
        quant=output_dir + "/salmon/counts_separate/{sample}/quant.sf"
    params:
        outdir=output_dir + "/salmon/counts_separate/{sample}",
        libtype="A",  # Automatic library type detection
        index=f"{genome_dir}/organisms/{organism}/SalmonIndex"
    log:
        "logs/salmon_quant.{sample}.log"
    benchmark:
        "benchmarks/salmon_quant.{sample}.txt"
    threads: threads
    conda: "envs/basic.yaml"
    shell:
        """
        salmon quant -i {params.index} -l {params.libtype} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} -o {params.outdir} --gcBias --validateMappings > {log} 2>&1
        """

rule salmon_promoter_counts:
    input:
        quant=output_dir + "/salmon/counts_separate/{sample}/quant.sf",
        promoter_rds = proactiv_rds
    output:
        counts=output_dir + "/salmon/counts_separate/{sample}_promoter_counts.rds",
    params:
        sample = lambda wildcards: wildcards.sample
    log:
        "logs/salmon_promoter_counts.{sample}.log"
    benchmark:
        "benchmarks/salmon_promoter_counts.{sample}.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/salmon_counts.R \
            {output_dir}/salmon/counts_separate \
            {input.promoter_rds} \
            {input.quant} \
            {params.sample} > {log} 2>&1
        """


rule salmon_promoter_merge:
    input:
        counts=expand(output_dir + "/salmon/counts_separate/{sample}_promoter_counts.rds", sample=samples),
        promoter_rds = proactiv_rds
    output:
        combined=output_dir + "/salmon/counts_merged/merged_promoter_counts.rds"
    params:
        samples=samples,
        condition = condition_str,
        condition_compare = condition_compare_str,
        newnames=samples,
    log:
          "logs/salmon_promoter_merge.log"    
    shell:
        """
        Rscript {workflow.basedir}/../scripts/merge_promoter_counts_copy.R \
            {output_dir}/salmon/counts_merged \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            "{params.condition_compare}"\
            {input.counts}  > {log} 2>&1
        """

rule salmon_promoter_merge_samplesheet:
    input:
        counts=expand(output_dir + "/salmon/counts_separate/{sample}_promoter_counts.rds", sample=diff_samples),
        promoter_rds = proactiv_rds
    output:
        combined = salmon_merge_out
    params:
        out_dir = salmon_merge_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        newnames = diff_samples_str,
    log:
          "logs/salmon_promoter_merge_samplesheet.log"  
    shell:
        """
        mkdir -p "{params.out_dir}"
        Rscript {workflow.basedir}/../scripts/merge_promoter_counts.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            "{params.condition_compare}"\
            {input.counts}  > {log} 2>&1
        """

rule salmon_promoter_classification_total:
    output:
        rowData = output_dir + "/salmon/promoter_classification_total/Summary_classified_rowData.rds",
        SE=output_dir + "/salmon/promoter_classification_total/Promoter_activity_SE.rds"
    input:
        combined=output_dir + "/salmon/counts_merged/merged_promoter_counts.rds",
        promoter_rds=proactiv_rds
    params:
        samples=" ".join(samples),
        newnames=" ".join(samples),
        condition = condition_str,
        condition_compare = condition_compare_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        batch=batch_condition,
        norm_method=config.get("norm_method", "deseq2")

    log:
        "logs/salmon_promoter_classify_total.log"
    benchmark:
        "benchmarks/salmon_promoter_classify_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/promoter_classify_total.R \
            {output_dir}/salmon/promoter_classification_total \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.combined} \
            "{params.condition_compare}" \
            "{params.fit_script}" \
            "{params.batch}" \
            "{params.norm_method}" > {log} 2>&1
        """


rule salmon_promoter_classification:
    input:
        combined=salmon_merge_out,
        promoter_rds=proactiv_rds
    output:
        se       = salmon_class_outputs["se"],
        rowData  = salmon_class_outputs["rowData"],
        raw_count= salmon_class_outputs["raw_count"]
    params:
        out_dir = salmon_class_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        newnames = diff_samples_str,
        batch_unsorted = diff_batch_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        test_condition = test_condition,
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/salmon_promoter_classify.log"
    benchmark:
        "benchmarks/salmon_promoter_classify.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_classify.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.combined} \
            "{params.condition_compare}"\
            "{params.fit_script}" \
            "{params.batch_unsorted}" \
            {params.norm_method} > {log} 2>&1
        """

rule salmon_differential:
    input:
        SE       = salmon_class_outputs["se"],
        promoter_rds = proactiv_rds,
        raw_count= salmon_class_outputs["raw_count"]
    output:
        deg = salmon_diff_out["DEG"],
        altpromoter    = salmon_diff_out["AP"]
    params:
        out_dir = salmon_diff_dir,
        samples       = diff_samples_str,
        condition = diff_conditions_str,
        baseline      = control_condition,
        reference    = test_condition,
        min_promoter_fold_change = min_pFC,
        max_gene_fold_change = max_gFC,
        lfcshrink = lfcshrink,
        batch = diff_batch_str,
    log:
        "logs/salmon_promoter_differential.log"
    benchmark:
        "benchmarks/salmon_promoter_differential.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_differential.R \
           "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.raw_count} \
            {input.SE} \
            {params.baseline} \
            {params.reference} \
            {params.min_promoter_fold_change} \
            {params.max_gene_fold_change} {params.lfcshrink} {params.batch} > {log} 2>&1 
        """


rule salmon_plots:
    input:
        SE       = salmon_class_outputs["se"],
        rowData  = salmon_class_outputs["rowData"],
        promoter_rds = proactiv_rds
    output:
        plots = salmon_plot_outputs
    params:
        out_dir = salmon_plots_dir,
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
    log:
        "logs/salmon_plots.log"
    benchmark:
        "benchmarks/salmon_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/plots.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """


salmon_overall_plot_outputs = [
    f"{output_dir}/salmon/plots_total/overall/overall_promoter_activity_{plot}.pdf"
    for plot in [
        "category_percentage", "category_comparison", "position_category",
        "geneexpression_correlation", "category_percentage_genewise",
        "single_multiple_category",  # Plot 1–6
        "number_hist_all", "number_hist_without1"  # Plot 7–9
    ]
]

rule salmon_overall_plots:
    input:
        rowData = output_dir + "/salmon/promoter_classification_total/Summary_classified_rowData.rds",
        promoter_rds = proactiv_rds,
        SE=output_dir + "/salmon/promoter_classification_total/Promoter_activity_SE.rds"
    output:
        plots = salmon_overall_plot_outputs
    params:
        condition = condition_str,
        condition_compare = condition_compare_str,
    log:
        "logs/salmon_plots_total.log"
    benchmark:
        "benchmarks/salmon_plots_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p {output_dir}/salmon/plots_total

        Rscript {workflow.basedir}/../scripts/plots_total.R \
            {output_dir}/salmon/plots_total \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """

rule method_comparison_plots:
    input:
        CAGE_gene = output_dir + "/../CAGE_HEART/salmon/compare/DESeq2_promoter_allResults_keepall.rds",
        CAGE_promoter = output_dir + "/../CAGE_HEART/salmon/compare/altPromoter_relativeUsage_delta10pct_allSamples1pct_fullAligned.rds",
        CAGE_minor = output_dir + "/../CAGE_HEART/salmon/plots_condition_wise/salmon_minor_promoter_class_list.rds",
        CAGE_major = output_dir + "/../CAGE_HEART/salmon/plots_condition_wise/salmon_major_promoter_class_list.rds",
        salmon_gene = output_dir + "/salmon/compare/DESeq2_promoter_allResults_keepall.rds",
        salmon_promoter = output_dir + "/salmon/compare/altPromoter_relativeUsage_delta10pct_allSamples1pct_fullAligned.rds",
        salmon_minor = output_dir + "/salmon/plots_condition_wise/salmon_minor_promoter_class_list.rds",
        salmon_major = output_dir + "/salmon/plots_condition_wise/salmon_major_promoter_class_list.rds",
        dexseq_gene = output_dir + "/dexseq/compare/DESeq2_promoter_allResults_keepall.rds",
        dexseq_promoter = output_dir + "/dexseq/compare/DESeq2_promoter_allResults_keepall.rds",
        dexseq_minor = output_dir + "/dexseq/plots_condition_wise/dexseq_minor_promoter_class_list.rds",
        dexseq_major = output_dir + "/dexseq/plots_condition_wise/dexseq_major_promoter_class_list.rds",
        proactiv_gene = output_dir + "/proactiv/compare/DESeq2_gene_allResults_keepall.rds",
        proactiv_promoter = output_dir + "/proactiv/compare/altPromoter_relativeUsage_delta10pct_allSamples1pct_fullAligned.rds",
        proactiv_category = output_dir + "/proactiv/quantify/rowData.rds",
    output:  
        major_venn = output_dir + "/comparison3/major_promoter_venn_combined.pdf",
        minor_venn = output_dir + "/comparison3/minor_promoter_venn_combined.pdf",
    #    promoter_up_FC = output_dir + "/comparison/log2FC_IQR_pval_up.pdf",
    #    promoter_down_FC = output_dir + "/comparison/log2FC_IQR_pval_down.pdf",
    #    gene_up_FC = output_dir + "/comparison/log2FC_promoter_IQR_FDR_up.pdf",
    #    gene_down_FC = output_dir + "/comparison/log2FC_promoter_IQR_FDR_down.pdf",
    #    promoter_reg_count = output_dir + "/comparison/promoter_count_pval.pdf",
    #    gene_reg_count = output_dir + "/comparison/promoter_count_FDR.pdf",
    params:
    log:
        "logs/method_comparison_plots.log"
    benchmark:
        "benchmarks/method_comparison_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/method_comparison_DEG_only.R \
            {output_dir}/comparison3 \
            {input.CAGE_gene} \
            {input.CAGE_promoter} \
            {input.CAGE_minor} \
            {input.CAGE_major} \
            {input.salmon_gene} \
            {input.salmon_promoter} \
            {input.salmon_minor} \
            {input.salmon_major} \
            {input.dexseq_gene} \
            {input.dexseq_promoter} \
            {input.dexseq_minor} \
            {input.dexseq_major} \
            {input.proactiv_gene} \
            {input.proactiv_promoter} \
            {input.proactiv_category} > {log} 2>&1
        """

#############################################
# CAGE Rules
#############################################
rule cage_featurecounts:
    input:
        #bam="CAGE_BAM/{sample}.bam",
        bam=output_dir + "/STAR/{sample}/{sample}.sorted.bam",
        gtf="/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/cage_promoter_counts.gtf"
    output:
        counts=output_dir + "/cage/counts_separate/{sample}_promoter_counts.txt",
        summary=output_dir + "/cage/counts_separate/{sample}_promoter_counts.txt.summary"
    log:
        "logs/counts_separate/{sample}.log"
    params:
        strand=1
    threads: threads
    conda: "envs/featurecounts.yaml"
    shell:
        """
        featureCounts \
          -a {input.gtf} \
          -o {output.counts} \
          -t promoter \
          -g promoter_id \
          -s {params.strand} -M -O \
          {input.bam} > {log} 2>&1
        """

rule cage_counts:
    input:
        promoter_rds = proactiv_rds,
        counts_txt = expand(output_dir + "/cage/counts_separate/{sample}_promoter_counts.txt", sample=samples)
    output:
        counts = output_dir + "/cage/counts_merged/cage_counts.rds"
    params:
        fc_dir = output_dir + "/cage/counts_separate",
        samples = " ".join(samples),
        newnames = " ".join(samples),
        condition = condition_str,
        condition_compare = condition_compare_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        batch = batch_condition
    log:
        "logs/cage_counts.log"
    benchmark:
        "benchmarks/cage_counts.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/CAGE_count.R \
            {output_dir}/cage/counts_merged \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {params.fc_dir} > {log} 2>&1
        """

rule cage_counts_diff:
    input:
        promoter_rds = proactiv_rds,
        counts_txt = expand(output_dir + "/cage/counts_separate/{sample}_promoter_counts.txt", sample=diff_samples)
    output:
        counts = cage_merge_out
    params:
        out_dir = cage_merge_dir,
        fc_dir = output_dir + "/cage/counts_separate",
        samples = " ".join(diff_samples),
        newnames = " ".join(diff_samples),
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        batch = batch_condition
    log:
        "logs/cage_counts.log"
    benchmark:
        "benchmarks/cage_counts.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir  -p "{params.out_dir}"
        Rscript {workflow.basedir}/../scripts/CAGE_count.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {params.fc_dir} > {log} 2>&1
        """

rule cage_promoter_classification_total:
    output:
        rowData = output_dir + "/cage/promoter_classification_total/Summary_classified_rowData.rds",
        SE=output_dir + "/cage/promoter_classification_total/Promoter_activity_SE.rds",
    input:
        counts=output_dir + "/cage/counts_merged/cage_counts.rds",
        promoter_rds=proactiv_rds
    params:
        fc_dir=output_dir + "/cage/featureCounts",
        samples=" ".join(samples),
        newnames=" ".join(samples),
        condition = condition_str,
        condition_compare = condition_compare_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        batch=batch_condition,
        norm_method=config.get("norm_method", "deseq2")

    log:
        "logs/cage_promoter_classification_total.log"
    benchmark:
        "benchmarks/cage_promoter_classification_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/promoter_classify_total.R \
            {output_dir}/cage/promoter_classification_total \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.counts} \
            "{params.condition_compare}" \
            "{params.fit_script}" \
            "{params.batch}" \
            "{params.norm_method}" > {log} 2>&1
        """

rule cage_promoter_classification:
    output:
        se       = cage_class_outputs["se"],
        rowData  = cage_class_outputs["rowData"],
        raw_count= cage_class_outputs["raw_count"]
    input:
        counts=cage_merge_out,
        promoter_rds=proactiv_rds
    params:
        out_dir = cage_class_dir,
        fc_dir=output_dir + "/cage/featureCounts",
        samples = " ".join(diff_samples),
        newnames = " ".join(diff_samples),
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        fit_script = os.path.join(SCRIPTS_DIR, "sanity_check.R"),
        batch=diff_batch_str,
        norm_method=config.get("norm_method", "deseq2")

    log:
        "logs/cage_promoter_classification.log"
    benchmark:
        "benchmarks/cage_promoter_classification.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_classify_cage.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.counts} \
            "{params.condition_compare}" \
            "{params.fit_script}" \
            "{params.batch}" \
            "{params.norm_method}" > {log} 2>&1
        """


rule cage_differential:
    input:
        SE       = cage_class_outputs["se"],
        promoter_rds  = proactiv_rds,
        raw_count= cage_class_outputs["raw_count"]
    output:
        deg = cage_diff_out["DEG"],
        altpromoter    = cage_diff_out["AP"]
    params:
        out_dir = cage_diff_dir,
        samples = " ".join(diff_samples),
        newnames = " ".join(diff_samples),
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
        batch = diff_batch_str,
        baseline      = control_condition,
        reference    = test_condition,
        min_promoter_fold_change = min_pFC,
        max_gene_fold_change = max_gFC,
        lfcshrink = lfcshrink,
    log:
        "logs/cage_promoter_differential.log"
    benchmark:
        "benchmarks/cage_promoter_differential.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/promoter_differential.R \
           "{params.out_dir}" \
            {input.promoter_rds} \
            "{params.samples}" \
            "{params.condition}" \
            {input.raw_count} \
            {input.SE} \
            {params.baseline} \
            {params.reference} \
            {params.min_promoter_fold_change} \
            {params.max_gene_fold_change} {params.lfcshrink} {params.batch} > {log} 2>&1 
        """

rule cage_plots_condition_wise:
    input:
        SE       = cage_class_outputs["se"],
        rowData  = cage_class_outputs["rowData"],
        promoter_rds  = proactiv_rds
    output:
        plots = cage_plot_outputs
    params:
        out_dir = cage_plots_dir,
        condition = diff_conditions_str,
        condition_compare = diff_condition_compare_str,
    log:
        "logs/cage_plots.log"
    benchmark:
        "benchmarks/cage_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p "{params.out_dir}"

        Rscript {workflow.basedir}/../scripts/plots.R \
            "{params.out_dir}" \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """

cage_overall_plot_outputs = [
    f"{output_dir}/cage/plots_total/overall/overall_promoter_activity_{plot}.pdf"
    for plot in [
        "category_percentage", "category_comparison", "position_category",
        "geneexpression_correlation", "category_percentage_genewise",
        "single_multiple_category",  # Plot 1–6
        "number_hist_all", "number_hist_without1"  # Plot 7–9
    ]
]

rule cage_overall_plots:
    input:
        rowData = output_dir + "/cage/promoter_classification_total/Summary_classified_rowData.rds",
        promoter_rds = proactiv_rds,
        SE=output_dir + "/cage/promoter_classification_total/Promoter_activity_SE.rds"
    output:
        plots = cage_overall_plot_outputs
    params:
        condition = condition_str,
        condition_compare = condition_compare_str,
    log:
        "logs/cage_plots_total.log"
    benchmark:
        "benchmarks/cage_plots_total.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p {output_dir}/cage/plots_total

        Rscript {workflow.basedir}/../scripts/plots_total.R \
            {output_dir}/cage/plots_total \
            {input.promoter_rds} \
            {input.rowData} \
            {input.SE} \
            "{params.condition_compare}" \
            "{params.condition}" > {log} 2>&1
        """

#############################################
# MultiQC Rule
#############################################
rule multiqc:
    input:
        *(salmon_overall_plot_outputs if do_salmon else []),
        *(proactiv_overall_plot_outputs if do_proactiv else []),
        *(dexseq_overall_plot_outputs if do_dexseq else []),
        *(cage_overall_plot_outputs if do_cage else []),
    output:
        report=output_dir + "/multiqc/multiqc_report.html"
    params:
        outdir=output_dir + "/multiqc",
        logdirs=[output_dir],
    log:
        "logs/multiqc.log"
    benchmark:
        "benchmarks/multiqc.txt"
    conda: "envs/featurecounts.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {params.logdirs} -o {params.outdir} --filename multiqc_report --force > {log} 2>&1
        """