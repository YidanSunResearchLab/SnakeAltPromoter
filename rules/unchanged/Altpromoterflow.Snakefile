import os

#############################################
# Configuration
#############################################
pipeline = config["pipeline"]
input_dir = config["input_dir"]
output_dir = config["output_dir"]
genome_dir = config["genome_dir"]
organism = config["organism"]
samples = config["samples"]  # Dictionary: {sample: {R1: path, R2: path}}
reads = ["R1", "R2"]  # Paired-end reads
downsample_size = config.get("downsample_size", 0)
threads = config.get("threads", 16)
trimmer_options = config.get("trimmer_options", "")
star_options = config.get("star_options", "")
samples_comparison = config["samples_comparison"]
samples_comparison = samples_comparison.split()

# Paths from genomesetup.smk
star_index = f"{genome_dir}/organisms/{organism}/STARIndex"
genes_gtf = f"{genome_dir}/organisms/{organism}/Annotation.gtf"
genes_bed = f"{genome_dir}/organisms/{organism}/Annotation/genes.bed"
dexseq_gff = f"{genome_dir}/organisms/{organism}/Annotation/DEXSeq_flattened_exons.gtf"
proactiv_rds = f"{genome_dir}/organisms/{organism}/Annotation/proActiv_promoter_annotation.rds"


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
rule all:
    input:
        # FastQC reports
        expand(output_dir + "/fastqc/{sample}_{read}_fastqc.html", sample=samples.keys(), read=reads),
        # Trimmed FASTQs
        expand(output_dir + "/fastqs/trimmed/{sample}_{read}.fastq.gz", sample=samples.keys(), read=reads),
        # STAR BAMs
        expand(output_dir + "/bam/{sample}.sorted.bam", sample=samples.keys()),
        # proActiv outputs
        #output_dir + "/proactiv/merged/normalized_promoter_counts.rds",
        #expand(output_dir + "/proactiv/plots/promoter_activity_{plot}.pdf", 
        #       plot=["category_percentage", "category_comparison", "position_category", 
        #             "geneexpression_correlation", "category_percentage_genewise", 
        #             "single_multiple_category", "number_hist_all", "number_hist_without1"]),
        # DEXSeq outputs
        #expand(output_dir + "/bam/{sample}.strand.txt", sample=samples.keys()),
        #expand(output_dir + "/dexseq/{sample}_counts.txt", sample=samples.keys()),
        #expand(output_dir + "/dexseq/{sample}_promoter_counts.rds", sample=samples.keys()),
        #output_dir + "/dexseq/normalized_promoter_counts.rds",
        #expand(output_dir + "/dexseq/plots/promoter_activity_{plot}.pdf", 
        #       plot=["category_percentage", "category_comparison", "position_category", 
        #             "geneexpression_correlation", "category_percentage_genewise", 
        #            "single_multiple_category", "number_hist_all", "number_hist_without1"])
        # Salmon outputs
        #expand(output_dir + "/salmon/{sample}/quant.sf", sample=samples.keys()),
        #output_dir + "/salmon/gene_counts.rds",
        #expand(output_dir + "/salmon/counts/{sample}_promoter_counts.rds", sample=samples.keys()),
        #output_dir + "/salmon/merged/normalized_promoter_counts.rds",
        #expand(output_dir + "/salmon/plots/promoter_activity_{plot}.pdf", 
        #       plot=["category_percentage", "category_comparison", "position_category", 
        #             "geneexpression_correlation", "category_percentage_genewise", 
        #             "single_multiple_category", "number_hist_all", "number_hist_without1"])

#############################################
# Preprocessing Rules
#############################################
rule link_fastqs:
    input:
        fastqs=input_dir+"/fastqs/raw/{sample}_{read}.fastq.gz"
    output:
        fastq=output_dir + "/fastqs/raw/{sample}_{read}.fastq.gz"
    log:
        "logs/link_fastqs.{sample}_{read}.log"
    benchmark:
        "benchmarks/link_fastqs.{sample}_{read}.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        ln -sf {input.fastqs} {output.fastq} 2> {log}
        """

rule fastqc:
    input:
        fastq=lambda wildcards: (
            output_dir + "/fastqs/downsampled/" + wildcards.sample + "_" + wildcards.read + ".fastq.gz"
            if downsample_size > 0
            else output_dir + "/fastqs/raw/" + wildcards.sample + "_" + wildcards.read + ".fastq.gz"
        )
    output:
        html=output_dir + "/fastqc/{sample}_{read}_fastqc.html",
        zip=output_dir + "/fastqc/{sample}_{read}_fastqc.zip"
    log:
        "logs/fastqc.{sample}_{read}.log"
    benchmark:
        "benchmarks/fastqc.{sample}_{read}.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        fastqc -o {output_dir}/fastqc {input.fastq} > {log} 2>&1
        """

rule downsample:
    input:
        fastq=output_dir + "/fastqs/raw/{sample}_{read}.fastq.gz"
    output:
        fastq=output_dir + "/fastqs/downsampled/{sample}_{read}.fastq.gz"
    params:
        size=downsample_size
    log:
        "logs/downsample.{sample}_{read}.log"
    benchmark:
        "benchmarks/downsample.{sample}_{read}.txt"
    threads: threads
    conda: "envs/basic.yaml"
    shell:
        """
        if [ {params.size} -gt 0 ]; then
            seqtk sample -s100 {input.fastq} {params.size} | gzip > {output.fastq} 2> {log}
        else
            ln -sf ../raw/{wildcards.sample}_{wildcards.read}.fastq.gz {output.fastq} 2> {log}
        fi
        """

rule trim_galore:
    input:
        r1=output_dir + "/fastqs/downsampled/{sample}_R1.fastq.gz",
        r2=output_dir + "/fastqs/downsampled/{sample}_R2.fastq.gz"
    output:
        r1=output_dir + "/fastqs/trimmed/{sample}_R1.fastq.gz",
        r2=output_dir + "/fastqs/trimmed/{sample}_R2.fastq.gz",
        report_r1=output_dir + "/fastqs/trimmed/{sample}_R1_fastqc.html",
        report_r2=output_dir + "/fastqs/trimmed/{sample}_R2_fastqc.html"
    params:
        opts=trimmer_options,
        outdir=output_dir + "/fastqs/trimmed"
    log:
        "logs/trim_galore.{sample}.log"
    benchmark:
        "benchmarks/trim_galore.{sample}.txt"
    threads: 4
    conda: "envs/basic.yaml"
    shell:
        """
        trim_galore --paired --stringency 3 {params.opts} \
            --output_dir {params.outdir} --gzip \
            --fastqc --fastqc_args "-o {params.outdir}" \
            {input.r1} {input.r2} > {log} 2>&1
        mv {params.outdir}/{wildcards.sample}_R1_val_1.fq.gz {output.r1}
        mv {params.outdir}/{wildcards.sample}_R2_val_2.fq.gz {output.r2}
        mv {params.outdir}/{wildcards.sample}_R1_val_1_fastqc.html {output.report_r1}
        mv {params.outdir}/{wildcards.sample}_R2_val_2_fastqc.html {output.report_r2}
        """

rule star:
    input:
        r1=output_dir + "/fastqs/trimmed/{sample}_R1.fastq.gz",
        r2=output_dir + "/fastqs/trimmed/{sample}_R2.fastq.gz",
        index=star_index + "/SAindex",
        gtf=genes_gtf
    output:
        bam=temp(output_dir + "/bam/{sample}.sorted.bam"),
        sj=output_dir + "/bam/{sample}/{sample}.SJ.out.tab"
    params:
        opts=star_options,
        prefix=output_dir + "/bam/{sample}/{sample}.",
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
        STAR --runThreadN {threads} \
            {params.opts} \
            --sjdbOverhang 100 \
            --outSAMunmapped Within \
            --outSAMtype BAM Unsorted \
            --outStd BAM_Unsorted \
            --sjdbGTFfile {input.gtf} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand 'gunzip -c' \
            --outFileNamePrefix {params.prefix} \
        | samtools sort -m {params.samsort_memory} -T {params.sample_dir}/{wildcards.sample} -@ {params.samtools_threads} -O bam -o {output.bam} - > {log} 2>&1
        """

#############################################
# proActiv Rules
#############################################
rule proactiv_counts:
    input:
        sj=output_dir + "/bam/{sample}/{sample}.SJ.out.tab",
        promoter_rds=proactiv_rds
    output:
        counts=output_dir + "/proactiv/counts/{sample}_promoter_counts.rds"
    params:
        sample="{sample}"
    log:
        "logs/proactiv_counts.{sample}.log"
    benchmark:
        "benchmarks/proactiv_counts.{sample}.txt"
    threads: threads
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/proactiv_counts.R \
            {output_dir}/proactiv/counts \
            {input.promoter_rds} \
            {input.sj} \
            {params.sample} > {log} 2>&1
        """

rule proactiv_merge:
    input:
        counts=expand(output_dir + "/proactiv/counts/{sample}_promoter_counts.rds", sample=samples.keys()),
        promoter_rds=proactiv_rds,
        sizefactor=output_dir + "/bam/dexseq_counts.txt"
    output:
        norm_counts=output_dir + "/proactiv/merged/normalized_promoter_counts.rds",
        abs_activity=output_dir + "/proactiv/merged/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/proactiv/merged/gene_expression.rds",
        rel_activity=output_dir + "/proactiv/merged/relative_promoter_activity.rds",
        abs_mean=output_dir + "/proactiv/merged/absolute_promoter_activity_mean.rds",
        rel_mean=output_dir + "/proactiv/merged/relative_promoter_activity_mean.rds",
        gene_mean=output_dir + "/proactiv/merged/gene_expression_mean.rds",
        category=output_dir + "/proactiv/merged/absolute_promoter_activity_category.rds"
    params:
        samples=list(samples.keys()),
        conditions=[config["conditions"].get(s, s) for s in samples.keys()],
        newnames=[config["sample_info"].get(s, s) for s in samples.keys()],
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/proactiv_merge.log"
    benchmark:
        "benchmarks/proactiv_merge.txt"
    threads: threads
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/proactiv_merge.R \
            {output_dir}/proactiv/merged \
            {input.promoter_rds} \
            {input.sizefactor} \
            {params.norm_method} \
            "{params.samples}" \
            "{params.conditions}" \
            "{params.newnames}" \
            {input.counts} > {log} 2>&1
        """

rule proactiv_plots:
    input:
        abs_activity=output_dir + "/proactiv/merged/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/proactiv/merged/gene_expression.rds",
        category=output_dir + "/proactiv/merged/absolute_promoter_activity_category.rds",
        promoter_rds=proactiv_rds
    output:
        expand(output_dir + "/proactiv/plots/promoter_activity_{plot}.pdf", 
               plot=["category_percentage", "category_comparison", "position_category", 
                     "geneexpression_correlation", "category_percentage_genewise", 
                     "single_multiple_category", "number_hist_all", "number_hist_without1"])
    params:
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/proactiv_plots.log"
    benchmark:
        "benchmarks/proactiv_plots.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/proactiv_plots.R \
            {output_dir}/proactiv/plots \
            {input.promoter_rds} \
            {input.abs_activity} \
            {input.gene_expr} \
            {input.category} \
            {params.norm_method} > {log} 2>&1
        """

#############################################
# DEXSeq Rules
#############################################
rule infer_strand:
    input:
        bam=output_dir + "/bam/{sample}.sorted.bam",
        bed=genes_bed
    output:
        strand=output_dir + "/bam/{sample}.strand.txt"
    log:
        "logs/infer_strand.{sample}.log"
    benchmark:
        "benchmarks/infer_strand.{sample}.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        infer_experiment.py -r {input.bed} -i {input.bam} > {output.strand} 2> {log}
        """

rule dexseq_featurecounts:
    input:
        bam=output_dir + "/bam/{sample}.sorted.bam",
        strand=output_dir + "/bam/{sample}.strand.txt",
        gtf=dexseq_gff
    output:
        counts=output_dir + "/dexseq/{sample}_counts.txt"
    params:
        threads=16
    log:
        "logs/dexseq_featurecounts.{sample}.log"
    benchmark:
        "benchmarks/dexseq_featurecounts.{sample}.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        strand=$(python {workflow.basedir}/../scripts/get_strand.py {input.strand})
        featureCounts -f -O -s $strand -p -T {params.threads} -F GTF \
            -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
        """

rule dexseq_counts:
    input:
        counts=output_dir + "/dexseq/{sample}_counts.txt",
        promoter_rds=proactiv_rds
    output:
        counts=output_dir + "/dexseq/{sample}_promoter_counts.rds"
    params:
        sample="{sample}",
        newname=lambda wildcards: config["sample_info"].get(wildcards.sample, wildcards.sample)
    log:
        "logs/dexseq_counts.{sample}.log"
    benchmark:
        "benchmarks/dexseq_counts.{sample}.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/dexseq_counts.R \
            {output_dir}/dexseq \
            {input.counts} \
            {input.promoter_rds} \
            {params.sample} \
            {params.newname} > {log} 2>&1
        """

rule dexseq_merge:
    input:
        counts=expand(output_dir + "/dexseq/{sample}_promoter_counts.rds", sample=samples.keys()),
        promoter_rds=proactiv_rds,
        feature_counts=output_dir + "/featureCounts/counts.tsv"
    output:
        norm_counts=output_dir + "/dexseq/normalized_promoter_counts.rds",
        abs_activity=output_dir + "/dexseq/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/dexseq/gene_expression.rds",
        rel_activity=output_dir + "/dexseq/relative_promoter_activity.rds",
        abs_mean=output_dir + "/dexseq/absolute_promoter_activity_mean.rds",
        rel_mean=output_dir + "/dexseq/relative_promoter_activity_mean.rds",
        gene_mean=output_dir + "/dexseq/gene_expression_mean.rds",
        category=output_dir + "/dexseq/absolute_promoter_activity_category.rds",
        sizefactor=output_dir + "/bam/dexseq_counts.txt"
    params:
        samples=list(samples.keys()),
        conditions=[config["conditions"].get(s, s) for s in samples.keys()],
        newnames=[config["sample_info"].get(s, s) for s in samples.keys()],
        norm_method=config.get("norm_method", "deseq2"),
        comparison=samples_comparison
    log:
        "logs/dexseq_merge.log"
    benchmark:
        "benchmarks/dexseq_merge.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/dexseq_merge.R \
            {output_dir}/dexseq \
            {input.promoter_rds} \
            {input.feature_counts} \
            {params.norm_method} \
            "{params.samples}" \
            "{params.conditions}" \
            "{params.newnames}" \
            "{params.comparison}" \
            {input.counts} > {log} 2>&1
        """

rule dexseq_plots:
    input:
        abs_activity=output_dir + "/dexseq/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/dexseq/gene_expression.rds",
        category=output_dir + "/dexseq/absolute_promoter_activity_category.rds",
        promoter_rds=proactiv_rds
    output:
        expand(output_dir + "/dexseq/plots/promoter_activity_{plot}.pdf", 
               plot=["category_percentage", "category_comparison", "position_category", 
                     "geneexpression_correlation", "category_percentage_genewise", 
                     "single_multiple_category", "number_hist_all", "number_hist_without1"])
    params:
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/dexseq_plots.log"
    benchmark:
        "benchmarks/dexseq_plots.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/proactiv_plots.R \
            {output_dir}/dexseq/plots \
            {input.promoter_rds} \
            {input.abs_activity} \
            {input.gene_expr} \
            {input.category} \
            {params.norm_method} > {log} 2>&1
        """

#############################################
# Salmon Rules
#############################################
rule salmon_quant:
    input:
        r1=output_dir + "/fastqs/trimmed/{sample}_R1.fastq.gz",
        r2=output_dir + "/fastqs/trimmed/{sample}_R2.fastq.gz",
    output:
        quant=output_dir + "/salmon/{sample}/quant.sf"
    params:
        outdir=output_dir + "/salmon/{sample}",
        libtype="A",  # Automatic library type detection
        index=f"{genome_dir}/organisms/{organism}/salmon_index"
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

rule salmon_merge:
    input:
        quants=expand(output_dir + "/salmon/{sample}/quant.sf", sample=samples.keys()),
        tx2gene=tx2gene
    output:
        counts=output_dir + "/salmon/gene_counts.rds"
    params:
        samples=list(samples.keys()),
        newnames=[config["sample_info"].get(s, s) for s in samples.keys()]
    log:
        "logs/salmon_merge.log"
    benchmark:
        "benchmarks/salmon_merge.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/salmon_merge.R \
            {output_dir}/salmon \
            {input.tx2gene} \
            "{params.samples}" \
            "{params.newnames}" \
            {input.quants} > {log} 2>&1
        """

rule salmon_promoter_counts:
    input:
        quant=output_dir + "/salmon/{sample}/quant.sf",
        promoter_rds=proactiv_rds
    output:
        counts=output_dir + "/salmon/counts/{sample}_promoter_counts.rds"
    params:
        sample="{sample}",
        newname=lambda wildcards: config["sample_info"].get(wildcards.sample, wildcards.sample)
    log:
        "logs/salmon_promoter_counts.{sample}.log"
    benchmark:
        "benchmarks/salmon_promoter_counts.{sample}.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/salmon_promoter_counts.R \
            {output_dir}/salmon/counts \
            {input.promoter_rds} \
            {input.quant} \
            {params.sample} \
            {params.newname} > {log} 2>&1
        """

rule salmon_promoter_merge:
    input:
        counts=expand(output_dir + "/salmon/counts/{sample}_promoter_counts.rds", sample=samples.keys()),
        promoter_rds=proactiv_rds,
        feature_counts=output_dir + "/featureCounts/counts.tsv"
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
        samples=list(samples.keys()),
        conditions=[config["conditions"].get(s, s) for s in samples.keys()],
        newnames=[config["sample_info"].get(s, s) for s in samples.keys()],
        norm_method=config.get("norm_method", "deseq2"),
        comparison=samples_comparison
    log:
        "logs/salmon_promoter_merge.log"
    benchmark:
        "benchmarks/salmon_promoter_merge.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/salmon_promoter_merge.R \
            {output_dir}/salmon/merged \
            {input.promoter_rds} \
            {input.feature_counts} \
            {params.norm_method} \
            "{params.samples}" \
            "{params.conditions}" \
            "{params.newnames}" \
            "{params.comparison}" \
            {input.counts} > {log} 2>&1
        """

rule salmon_promoter_plots:
    input:
        abs_activity=output_dir + "/salmon/merged/absolute_promoter_activity.rds",
        gene_expr=output_dir + "/salmon/merged/gene_expression.rds",
        category=output_dir + "/salmon/merged/absolute_promoter_activity_category.rds",
        promoter_rds=proactiv_rds
    output:
        expand(output_dir + "/salmon/plots/promoter_activity_{plot}.pdf", 
               plot=["category_percentage", "category_comparison", "position_category", 
                     "geneexpression_correlation", "category_percentage_genewise", 
                     "single_multiple_category", "number_hist_all", "number_hist_without1"])
    params:
        norm_method=config.get("norm_method", "deseq2")
    log:
        "logs/salmon_promoter_plots.log"
    benchmark:
        "benchmarks/salmon_promoter_plots.txt"
    conda: "envs/basic.yaml"
    shell:
        """
        Rscript {workflow.basedir}/../scripts/proactiv_plots.R \
            {output_dir}/salmon/plots \
            {input.promoter_rds} \
            {input.abs_activity} \
            {input.gene_expr} \
            {input.category} \
            {params.norm_method} > {log} 2>&1
        """