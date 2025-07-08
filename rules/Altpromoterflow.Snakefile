import os
print("CONFIG = ", config)
print("METHOD = ", config.get("method", "not provided"))
#############################################
# Configuration
#############################################
pipeline = config["pipeline"]
input_dir = config["input_dir"]
output_dir = config["output_dir"]
genome_dir = config["genome_dir"]
organism = config["organism"]
samples = config["samples"]  # Not dictionary but list # NOT Dictionary: {sample: {R1: path, R2: path}}

# Separate R1 and R2 so they are not read as one read
if isinstance(config["reads"], str):
    reads = config["reads"].replace(",", " ").split()
else:
    reads = config["reads"]
print("READS passed to Snakemake:", reads)

is_paired = "R2" in reads

#Differ rules for single-ended or paired reads due to differing outputs or different modes used
if is_paired:
    ruleorder: dexseq_featurecounts_paired > dexseq_featurecounts_single
    ruleorder: trim_galore_paired > trim_galore_single
    ruleorder: star_paired > star_single
    ruleorder: salmon_quant_paired > salmon_quant_single
else:
    ruleorder: dexseq_featurecounts_single > dexseq_featurecounts_paired
    ruleorder: trim_galore_single > trim_galore_paired
    ruleorder: star_single > star_paired
    ruleorder: salmon_quant_single > salmon_quant_paired


downsample_size = config.get("downsample_size", 0)
threads = config.get("threads", 16)
trimmer_options = config.get("trimmer_options", "")
star_options = config.get("star_options", "")
reference_condition = config.get("reference_condition","")


# Choose method to continue analysis
method = config.get("method", "all").lower()
do_salmon = method in ["salmon", "all"]
do_proactiv = method in ["proactiv", "all"]
do_dexseq = method in ["dexseq", "all"]

import re

# Use unsorted samples in salmon and dexseq
sample_to_group_unsorted = {}
batch_unsorted = []

# "_" is omitted if one column is blank during naming to prevent mutiple consecutive "_" (e.g., "__"). 
# The batch we are using have more blanks for single-ended, thus the part where condition is at differs from paired reads.
for s in samples:
    parts = s.split("_")
    if is_paired:
        condition = parts[6]  # e.g., "Ctl" or "Zta"
    else:
        condition = parts[2]
    sample_to_group_unsorted[s] = condition
    batch_unsorted.append("batch1")

samples_comparison_unsorted = [sample_to_group_unsorted[s] for s in samples]
condition_unsorted = ",".join(samples_comparison_unsorted)
unsorted_comparison = ",".join(sorted(set(samples_comparison_unsorted)))
batch_unsorted_str = ",".join(batch_unsorted)
print("Condition_unsorted passed to Snakemake:", condition_unsorted)

# Need sorted_samples to match STAR file order for proactiv
sorted_samples = sorted(samples)

sample_to_group = {s: sample_to_group_unsorted[s] for s in sorted_samples}
samples_comparison = [sample_to_group[s] for s in sorted_samples]
condition = ",".join(samples_comparison)
batch = [batch_unsorted[samples.index(s)] for s in sorted_samples]
batch_str = ",".join(batch)

# Get unique cell lines in order of appearance
cell_line = list(dict.fromkeys(samples_comparison))
cell_line_str = " ".join(cell_line)

# Paths from genomesetup.smk
star_index = f"{genome_dir}/organisms/{organism}/STARIndex"
genes_gtf = f"{genome_dir}/organisms/{organism}/Annotation.gtf"
genes_bed = f"{genome_dir}/organisms/{organism}/Annotation/genes.bed"
dexseq_gff = f"{genome_dir}/organisms/{organism}/Annotation/DEXSeq_flattened_exons.gff"
proactiv_rds = f"{genome_dir}/organisms/{organism}/Annotation/proActiv_promoter_annotation.rds"
tx2gene=f"{genome_dir}/organisms/{organism}/Annotation/genes_t2g.tsv"


# Validate config
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
os.makedirs(output_dir + "/dexseq/plots", exist_ok=True)
os.makedirs(output_dir + "/salmon", exist_ok=True)
os.makedirs(output_dir + "/salmon/counts", exist_ok=True)
os.makedirs(output_dir + "/salmon/merged", exist_ok=True)
os.makedirs(output_dir + "/salmon/plots", exist_ok=True)
os.makedirs("logs", exist_ok=True)
os.makedirs("benchmarks", exist_ok=True)


#############################################
# Final Target
#############################################
input_all = []

# Salmon outputs
if do_salmon:
    input_all += expand(output_dir + "/salmon/{sample}/quant.sf", sample=samples)
    input_all += expand(output_dir + "/salmon/counts/{sample}_promoter_counts.rds", sample=samples)
    input_all.append(output_dir + "/salmon/merged/normalized_promoter_counts.rds")
    input_all += expand(output_dir + "/salmon/plots/promoter_activity_{plot}.pdf", 
                        plot=["number_hist_all", "number_hist_without1", "tsne_plot"])
    input_all += expand(output_dir + "/salmon/plots/{cell_line}/{cell_line}_promoter_activity_{plot}.pdf", 
                        cell_line=cell_line, 
                        plot=["category_comparison", "geneexpression_correlation", "position_category", 
                              "category_percentage_genewise", "category_percentage", "single_multiple_category"])

# proActiv outputs
if do_proactiv:
    input_all.append(output_dir + "/proactiv/quantify/proactiv_result.rds")
    input_all += expand(output_dir + "/proactiv/plots/promoter_activity_{plot}.pdf", 
                        plot=["number_hist_all", "number_hist_without1", "tsne_plot"])
    input_all += expand(output_dir + "/proactiv/plots/{cell_line}/{cell_line}_promoter_activity_{plot}.pdf", 
                        cell_line=cell_line, 
                        plot=["category_comparison", "geneexpression_correlation", "position_category", 
                              "category_percentage_genewise", "category_percentage", "single_multiple_category"])

# DEXSeq outputs
if do_dexseq:
    input_all += expand(output_dir + "/bam/{sample}/{sample}.strand.txt", sample=samples)
    input_all += expand(output_dir + "/dexseq/{sample}/{sample}_counts.txt", sample=samples)
    input_all += expand(output_dir + "/dexseq/{sample}/{sample}_promoter_counts.rds", sample=samples)
    input_all.append(output_dir + "/dexseq/merge/normalized_promoter_counts.rds")
    input_all.append(output_dir + "/dexseq/merge/goodness_of_fit_ks_pvalues.rds")
    input_all.append(output_dir + "/dexseq/merge/sanity_check_DESeq2_result.rds")
    input_all += expand(output_dir + "/dexseq/plots/promoter_activity_{plot}.pdf", 
                        plot=["number_hist_all", "number_hist_without1", "tsne_plot"])
    input_all += expand(output_dir + "/dexseq/plots/{cell_line}/{cell_line}_promoter_activity_{plot}.pdf", 
                        cell_line=cell_line, 
                        plot=["category_comparison", "geneexpression_correlation", "position_category", 
                              "category_percentage_genewise", "category_percentage", "single_multiple_category"])


rule all:
    input:
        # FastQs link
        expand(output_dir + "/fastqs/raw/{sample}/{sample}_{read}.fastq.gz", sample=samples, read=reads),
        # FastQC reports
        expand(output_dir + "/fastqc/{sample}/{sample}_{read}_fastqc.html", sample=samples, read=reads),
        expand(output_dir + "/fastqc/{sample}/{sample}_{read}_fastqc.zip", sample=samples, read=reads),
        # Trimmed FASTQs
        expand(output_dir + "/fastqs/trimmed/{sample}/{sample}_{read}.fastq.gz", sample=samples, read=reads),
        # STAR BAMs
        expand(output_dir + "/bam/{sample}/{sample}.sorted.bam", sample=samples),
        # Based on method choice
        input_all

#############################################
# Preprocessing Rules
#############################################
rule link_fastqs:
    input:
        fastqs=input_dir+"/{sample}_{read}.fastq.gz"
    output:
        fastq=output_dir + "/fastqs/raw/{sample}/{sample}_{read}.fastq.gz"
    log:
        output_dir + "/logs/link_fastqs.{sample}_{read}.log"
    benchmark:
        "benchmarks/link_fastqs.{sample}_{read}.txt"
    conda: "envs/altbasic.yaml"
    shell:
        """
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
    benchmark:
        "benchmarks/fastqc.{sample}_{read}.txt"
    conda: "envs/altbasic.yaml"
    shell:
        """
        mkdir -p {output_dir}/fastqc/{wildcards.sample}
        base=$(basename {input.fastq} .fastq.gz)
        fastqc -o {output_dir}/fastqc {input.fastq} > {log} 2>&1
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
    threads: 8
    conda: "envs/altbasic.yaml"
    shell:
        """
        trim_galore --stringency 3 {params.opts} \
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
    threads: 4
    conda: "envs/altbasic.yaml"
    shell:
        """
        mkdir -p $(dirname {output.r1}) $(dirname {output.r2})

        trim_galore --paired --stringency 3 {params.opts} \
            --output_dir {params.outdir} --gzip \
            --fastqc --fastqc_args "-o {params.outdir}" \
            {input.r1} {input.r2} > {log} 2>&1
        mv {params.outdir}/{wildcards.sample}_R1_val_1.fq.gz {output.r1}
        mv {params.outdir}/{wildcards.sample}_R2_val_2.fq.gz {output.r2}
        mv {params.outdir}/{wildcards.sample}_R1_val_1_fastqc.html {output.report_r1}
        mv {params.outdir}/{wildcards.sample}_R2_val_2_fastqc.html {output.report_r2}
        """

rule star_single:
    input:
        r1=output_dir + "/fastqs/trimmed/{sample}/{sample}_R1.fastq.gz",
        index=star_index,
        gtf=genes_gtf
    output:
        bam=temp(output_dir + "/bam/{sample}/{sample}.sorted.bam"),
        sj=output_dir + "/bam/junctions/{sample}.SJ.out.tab"
    params:
        opts=star_options,
        prefix=output_dir + "/bam/junctions/{sample}.",
        sample_dir=output_dir + "/bam/{sample}",
        samsort_memory="8G",
        samtools_threads=16
    log:
        "logs/star.{sample}.log"
    benchmark:
        "benchmarks/star.{sample}.txt"
    threads: min(20, threads)
    conda: "envs/basic.yaml"
    shell:
        """
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
            --outFileNamePrefix {params.prefix} \
        
         samtools sort -m {params.samsort_memory} -T {params.sample_dir}/{wildcards.sample} -@ {params.samtools_threads} -O bam -o {output.bam} {params.prefix}Aligned.out.bam > {log} 2>&1
        """

rule star_paired:
    input:
        r1=output_dir + "/fastqs/trimmed/{sample}/{sample}_R1.fastq.gz",
        r2=output_dir + "/fastqs/trimmed/{sample}/{sample}_R2.fastq.gz",
        index=star_index,
        gtf=genes_gtf
    output:
        bam=temp(output_dir + "/bam/{sample}/{sample}.sorted.bam"),
        sj=output_dir + "/bam/junctions/{sample}.SJ.out.tab"
    params:
        opts=star_options,
        prefix=output_dir + "/bam/junctions/{sample}.",
        sample_dir=output_dir + "/bam/{sample}",
        samsort_memory="8G",
        samtools_threads=16
    log:
        "logs/star.{sample}.log"
    benchmark:
        "benchmarks/star.{sample}.txt"
    threads: min(20, threads)
    conda: "envs/basic.yaml"
    shell:
        """
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
            --outFileNamePrefix {params.prefix} \
        
         samtools sort -m {params.samsort_memory} -T {params.sample_dir}/{wildcards.sample} -@ {params.samtools_threads} -O bam -o {output.bam} {params.prefix}Aligned.out.bam > {log} 2>&1
        """

#############################################
# proActiv Rules
#############################################

rule proactiv_quantify:
    output:
        quantified=output_dir + "/proactiv/quantify/proactiv_result.rds",
        rowData = output_dir + "/proactiv/quantify/rowData.rds",
        gene_expr = output_dir + "/proactiv/quantify/gene_expression.rds"
    input:
        promoter_rds=proactiv_rds,
        sj_files=expand(output_dir + "/bam/junctions/{sample}.SJ.out.tab", sample=samples)
    params:
        condition = condition,
        reference_condition = reference_condition,
        batch = batch_str,
        fit_script="../scripts/sanity_check.R",
    log:
        "logs/proactiv_quantify.log"
    benchmark:
        "benchmarks/proactiv_quantify.txt"
    threads: threads
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p {output_dir}/proactiv/quantify

        FILES=$(echo {input.sj_files} | tr ' ' '\n' | sort | paste -sd' ')
        echo "FILES=$FILES"
        echo "CONDITION={params.condition}"
        Rscript {workflow.basedir}/../scripts/proactiv_quantify.R \
            {output_dir}/proactiv/quantify \
            {input.promoter_rds} \
            "$FILES" \
            {params.condition} \
            {params.reference_condition}\
            "{params.batch}" \
            "{params.fit_script}" > {log} 2>&1
        """

rule proactiv_plots:
    input:
        rowData = output_dir + "/proactiv/quantify/rowData.rds",
        gene_expr = output_dir + "/proactiv/quantify/gene_expression.rds",
        promoter_rds = proactiv_rds,
        result = output_dir + "/proactiv/quantify/proactiv_result.rds"
    output:
        expand(output_dir + "/proactiv/plots/promoter_activity_{plot}.pdf",
            plot=["number_hist_all", "number_hist_without1", "tsne_plot"]) +
        [
            f"{output_dir}/proactiv/plots/{cl}/{cl}_promoter_activity_{plot}.pdf"
            for cl in cell_line
            for plot in ["geneexpression_correlation", "category_percentage_genewise", "position_category", "category_percentage", "single_multiple_category", "category_comparison"]
        ]
    params:
        cell_line = cell_line_str,
        condition = condition
    log:
        "logs/proactiv_plots.log"
    benchmark:
        "benchmarks/proactiv_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        mkdir -p {output_dir}/proactiv/plots

        Rscript {workflow.basedir}/../scripts/proactiv_plots.R \
            {output_dir}/proactiv/plots \
            {input.promoter_rds} \
            {input.rowData} \
            {input.result} \
            {input.gene_expr} \
            "{params.cell_line}" \
            "{params.condition}" > {log} 2>&1
        """

#############################################
# DEXSeq Rules
#############################################
rule infer_strand:
    input:
        bam=output_dir + "/bam/{sample}/{sample}.sorted.bam",
        bed=genes_bed
    output:
        strand=output_dir + "/bam/{sample}/{sample}.strand.txt"
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
        bam=output_dir + "/bam/{sample}/{sample}.sorted.bam",
        strand=output_dir + "/bam/{sample}/{sample}.strand.txt",
        gff=dexseq_gff
    output:
        counts=output_dir + "/dexseq/{sample}/{sample}_counts.txt"
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
        bam=output_dir + "/bam/{sample}/{sample}.sorted.bam",
        strand=output_dir + "/bam/{sample}/{sample}.strand.txt",
        gff=dexseq_gff
    output:
        counts=output_dir + "/dexseq/{sample}/{sample}_counts.txt"
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

rule dexseq_counts:
    output:
        counts=output_dir + "/dexseq/{sample}/{sample}_promoter_counts.rds"
    input:
        gff=dexseq_gff,
        counts=output_dir + "/dexseq/{sample}/{sample}_counts.txt",
        promoter_rds=proactiv_rds
    params:
        samples="{sample}",
        newnames="{sample}"
    log:
        "logs/dexseq_counts.{sample}.log"
    benchmark:
        "benchmarks/dexseq_counts.{sample}.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/dexseq_counts.R \
            {output_dir}/dexseq \
            {input.gff} \
            {input.counts} \
            {input.promoter_rds} \
            {params.samples} \
            {params.newnames} > {log} 2>&1
        """

rule dexseq_merge:
    input:
        counts=expand(output_dir + "/dexseq/{sample}/{sample}_promoter_counts.rds", sample=samples),
        promoter_rds=proactiv_rds,
    output:
        norm_counts=output_dir + "/dexseq/merge/normalized_promoter_counts.rds",
        abs_activity=output_dir + "/dexseq/merge/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/dexseq/merge/gene_expression.rds",
        rel_activity=output_dir + "/dexseq/merge/relative_promoter_activity.rds",
        abs_mean=output_dir + "/dexseq/merge/absolute_promoter_activity_mean.rds",
        rel_mean=output_dir + "/dexseq/merge/relative_promoter_activity_mean.rds",
        gene_mean=output_dir + "/dexseq/merge/gene_expression_mean.rds",
        category=output_dir + "/dexseq/merge/absolute_promoter_activity_category.rds",
        sizefactor=output_dir + "/bam/dexseq_counts.txt",
        goodness_of_fit=output_dir + "/dexseq/merge/goodness_of_fit_ks_pvalues.rds",
        sanity_check=output_dir + "/dexseq/merge/sanity_check_DESeq2_result.rds",
    params:
        samples=samples,
        condition_unsorted = condition_unsorted,
        samples_comparison = unsorted_comparison,
        newnames=samples,
        batch_unsorted = batch_unsorted_str,
        fit_script="../scripts/sanity_check.R",
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/dexseq_merge.log"
    benchmark:
        "benchmarks/dexseq_merge.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/dexseq_merge.R \
            {output_dir}/dexseq/merge \
            {input.promoter_rds} \
            {params.norm_method} \
            "{params.samples}" \
            "{params.condition_unsorted}" \
            "{params.newnames}" \
            "{params.samples_comparison}" \
            "{params.batch_unsorted}" \
            "{params.fit_script}" \
            {input.counts} > {log} 2>&1
        """

rule dexseq_plots:
    input:
        abs_activity=output_dir + "/dexseq/merge/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/dexseq/merge/gene_expression.rds",
        promoter_rds=proactiv_rds
    output:
        expand(output_dir + "/dexseq/plots/promoter_activity_{plot}.pdf",
            plot=["number_hist_all", "number_hist_without1", "tsne_plot"]) +
        [
            f"{output_dir}/dexseq/plots/{cl}/{cl}_promoter_activity_{plot}.pdf"
            for cl in cell_line
            for plot in ["geneexpression_correlation", "category_percentage_genewise", "position_category", "category_percentage", "single_multiple_category", "category_comparison"]
        ]
    params:
        cell_lines=" ".join(cell_line),
        condition_unsorted = condition_unsorted,
    log:
        "logs/dexseq_plots.log"
    benchmark:
        "benchmarks/dexseq_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/dexseq_plots.R \
            {output_dir}/dexseq/plots \
            {input.promoter_rds} \
            {input.abs_activity} \
            {input.gene_expr} \
            "{params.cell_lines}"\
            "{params.condition_unsorted}" > {log} 2>&1
        """

#############################################
# Salmon Rules
#############################################
rule salmon_quant_single:
    input:
        r1=output_dir + "/fastqs/trimmed/{sample}/{sample}_R1.fastq.gz"
    output:
        quant=output_dir + "/salmon/{sample}/quant.sf"
    params:
        outdir=output_dir + "/salmon/{sample}",
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
        r1=output_dir + "/fastqs/trimmed/{sample}/{sample}_R1.fastq.gz",
        r2=output_dir + "/fastqs/trimmed/{sample}/{sample}_R2.fastq.gz",
    output:
        quant=output_dir + "/salmon/{sample}/quant.sf"
    params:
        outdir=output_dir + "/salmon/{sample}",
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
        quant=output_dir + "/salmon/{sample}/quant.sf",
        promoter_rds = proactiv_rds
    output:
        counts=output_dir + "/salmon/counts/{sample}_promoter_counts.rds"
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
            {output_dir}/salmon/counts \
            {input.promoter_rds} \
            {input.quant} \
            {params.sample} > {log} 2>&1
        """

rule salmon_promoter_merge:
    input:
        counts=expand(output_dir + "/salmon/counts/{sample}_promoter_counts.rds", sample=samples),
        promoter_rds=proactiv_rds
        
    output:
        norm_counts=output_dir + "/salmon/merged/normalized_promoter_counts.rds",
        abs_activity=output_dir + "/salmon/merged/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/salmon/merged/gene_expression.rds",
        rel_activity=output_dir + "/salmon/merged/relative_promoter_activity.rds",
        abs_mean=output_dir + "/salmon/merged/absolute_promoter_activity_mean.rds",
        rel_mean=output_dir + "/salmon/merged/relative_promoter_activity_mean.rds",
        gene_mean=output_dir + "/salmon/merged/gene_expression_mean.rds",
        category=output_dir + "/salmon/merged/absolute_promoter_activity_category.rds"
    params:
        samples=samples,
        condition_unsorted = condition_unsorted,
        newnames=samples,
        batch_unsorted = batch_unsorted_str,
        fit_script="../scripts/sanity_check.R",
        norm_method=config.get("norm_method", "deseq2")
        #norm_method="edger"
    log:
        "logs/salmon_promoter_merge.log"
    benchmark:
        "benchmarks/salmon_promoter_merge.txt"
    conda: "envs/altbasicR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/salmon_merge.R \
            {output_dir}/salmon/merged \
            {input.promoter_rds} \
            {params.norm_method} \
            "{params.samples}" \
            "{params.condition_unsorted}" \
            "{params.newnames}" \
            "{params.batch_unsorted}" \
            "{params.fit_script}" \
            {input.counts} > {log} 2>&1
        """

rule salmon_promoter_plots:
    input:
        abs_activity=output_dir + "/salmon/merged/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/salmon/merged/gene_expression.rds",
        promoter_rds=proactiv_rds
    output:
        expand(output_dir + "/salmon/plots/promoter_activity_{plot}.pdf",
            plot=["number_hist_all", "number_hist_without1", "tsne_plot"]) +
        [
            f"{output_dir}/salmon/plots/{cl}/{cl}_promoter_activity_{plot}.pdf"
            for cl in cell_line
            for plot in ["geneexpression_correlation", "category_percentage_genewise", "position_category", "category_percentage", "single_multiple_category", "category_comparison"]
        ]
    params:
        cell_lines=" ".join(cell_line),
        condition_unsorted = condition_unsorted
    log:
        "logs/salmon_promoter_plots.log"
    benchmark:
        "benchmarks/salmon_promoter_plots.txt"
    conda: "envs/dexR.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/salmon_plots.R \
            {output_dir}/salmon/plots \
            {input.promoter_rds} \
            {input.abs_activity} \
            {input.gene_expr} \
            "{params.cell_lines}"\
            "{params.condition_unsorted}" > {log} 2>&1
        """
