### `docs/Altpromoterflow.md`

# Altpromoterflow: Promoter Activity and Differential Analysis

`Altpromoterflow` maps reads onto genome, quantifies promoter activity, and compares across conditions using Trim Galore, STAR, proActiv, Salmon, and DEXSeq scripts.


## Usage
```bash
Altpromoterflow [options]
```

### Required Arguments
- `-i/--input_dir`: Directory containing FASTQ gz files.
- `-o/--output_dir`: Output directory.
- `--organism`: Genome assembly, the same as Genomesetup step (e.g., `hg38`).
- `--genome_dir`: Directory from Genomesetup step (e.g., `/abs/path/to/genome`).

### Optional Arguments
- `--pipeline`: Analysis software to use. Options: Maps, Fithichip, Hichipper, Hicdcplus, All. Default: Maps.
- `--downsample_size`: Number of valid pairs to downsample to (e.g., 50000000). Set to 0 to disable (default: 0).
- `--restriction_enzyme`: Restriction enzyme used in the experiment. Options: mboi (default), hindiii, dpnii, bglii, ncoi, msei, hinfI, mnase, arima.
- `--bin_size`: Interaction bin size in base pairs (e.g., 10000). Default: 5000.
- `--binning_range`: Maximum distance for binning interactions in base pairs (e.g., 1000000). Default: 2000000.
- `--length_cutoff`: Minimum fragment length for long contacts in base pairs (e.g., 1000). Default: 1000.
- `--threads`: Number of CPU threads for parallel processing (e.g., 100). Default: 30.
- `--fdr`: Minimum FDR threshold for significant interactions (e.g., 0.05). Default: 0.01.
- `--macs2_peaks`: Optional path to a MACS2 peaks file (e.g., peaks.narrowPeak) from ChIP-seq data to refine interactions. Default: NA (none).
- `--fithichip_BiasType`: FitHiChIP bias correction type. Options: 1 (coverage-based, default), 2 (ICE-based).
- `--maps_count_cutoff`: Minimum read count for MAPS significant interactions (e.g., 5). Default: 5.
- `--maps_ratio_cutoff`: Minimum observed-to-expected ratio for MAPS (e.g., 2.0). Default: 2.0.
- `--maps_model`: Statistical model for MAPS regression. Options: pospoisson (default), negbinom.
- `--maps_sex_chroms`: Sex chromosomes to include in MAPS analysis. Options: NA (none), X (default), Y, XY.
- `--hicpro_params`: Optional HiC-Pro config file path. Default: NA (none).
- `--hichipper_params`: Optional space-separated Hichipper parameters (e.g., --read-length 75). Default: NA (none).
- `--hicdc_params`: Optional space-separated Hicdcplus parameters (e.g., --PeakFile peaks.bed). Default: NA (none).
- `--samples_comparison`: Space-separated sample names for differential analysis (e.g., SampleA {sample}). Default: none.
- `--reference_condition`: Sample name for reference condition in differential analysis (e.g., SampleA). Default: NA (none).
- `--method`: Method to be used for differential analysis in Maps pipeline. Options are proactiv, salmon, dexseq, all. Default: all (all three methods).


### Example
```bash
Altpromoterflow.py \
  -i /path/to/your/fastqs/ \
  -o /path/to/your/output/ \
  --genome_dir /path/to/your/genome/ \
  --organism hg38 \
  --reference_condition SampleA \
  --threads 30 \
  --bin_size 5000 \
  --downsample_size 50000000 \
  --fdr 0.01 \
  --maps_model negbinom \
  --maps_sex_chroms X
```

## Altpromoterflow Output
```
tree AltPromoterFlow/
AltPromoterFlow/
├── bam
│   ├── {sample}
│   │   ├── {sample}.sorted.bam
│   │   └── {sample}.strand.txt
├── benchmarks
│   ├── dexseq_counts.{sample}.txt
│   ├── dexseq_featurecounts.{sample}.txt
│   ├── dexseq_fit.txt
│   ├── dexseq_merge.txt
│   ├── dexseq_plots.txt
│   ├── downsample.{sample}_R1.txt
│   ├── downsample.{sample}_R2.txt
│   ├── fastqc.{sample}_R1.txt
│   ├── fastqc.{sample}_R2.txt
│   ├── infer_strand.{sample}.txt
│   ├── link_fastqs.{sample}_R1.txt
│   ├── link_fastqs.{sample}_R2.txt
│   ├── proactiv_plots_addition.txt
│   ├── proactiv_plots.txt
│   ├── proactiv_quantify.txt
│   ├── salmon_promoter_counts.{sample}.txt
│   ├── salmon_promoter_merge.txt
│   ├── salmon_promoter_plots.txt
│   ├── salmon_quant.{sample}.txt
│   ├── star.{sample}.txt
│   └── trim_galore.{sample}.txt
├── dexseq
│   ├── {sample}
│   │   ├── {sample}_counts.txt
│   │   ├── {sample}_counts.txt.summary
│   │   └── {sample}_promoter_counts.rds
│   ├── dexseq.exon.bed
│   ├── dexseq.exoninfo.promoter.correspond.rds
│   ├── dexseq.exoninfo.promoter.correspond.txt
│   ├── dexseq.exoninfo.rds
│   ├── dexseq.transcripts.rds
│   ├── merge
│   │   ├── absolute_promoter_activity_category.rds
│   │   ├── absolute_promoter_activity_mean.rds
│   │   ├── absolute_promoter_activity.rds
│   │   ├── dispersion_per_gene.rds
│   │   ├── gene_expression_mean.rds
│   │   ├── gene_expression.rds
│   │   ├── goodness_of_fit_ks_pvalues.rds
│   │   ├── NB_size_per_gene.rds
│   │   ├── normalized_promoter_counts.rds
│   │   ├── per_sample_pnbinom_pvalues_matrix.rds
│   │   ├── promoter_counts.clean.rds
│   │   ├── promoter_counts.rds
│   │   ├── promoter_fitted_mu.rds
│   │   ├── relative_promoter_activity_mean.rds
│   │   ├── relative_promoter_activity.rds
│   │   ├── sanity_check_DESeq2_result.rds
│   │   └── size_factors.rds
│   ├── plots
│   │   ├── {cell_line}
│   │   │   ├── {cell_line}_promoter_activity_category_comparison.pdf
│   │   │   ├── {cell_line}_promoter_activity_category_percentage_genewise.pdf
│   │   │   ├── {cell_line}_promoter_activity_category_percentage.pdf
│   │   │   ├── {cell_line}_promoter_activity_geneexpression_correlation.pdf
│   │   │   ├── {cell_line}_promoter_activity_position_category.pdf
│   │   │   └── {cell_line}_promoter_activity_single_multiple_category.pdf
│   │   ├── promoter_activity_number_hist_all.pdf
│   │   ├── promoter_activity_number_hist_without1.pdf
│   │   └── promoter_activity_tsne_plot.pdf
├── fastqc
│   └── {sample}
│       ├── {sample}_R1_fastqc.html
│       ├── {sample}_R1_fastqc.zip
│       ├── {sample}_R2_fastqc.html
│       └── {sample}_R2_fastqc.zip
├── fastqs
│   ├── downsampled
│   │   └── {sample}
│   │   │   ├── {sample}_R1.fastq.gz -> 
│   │   │   └── {sample}_R2.fastq.gz -> 
│   ├── raw
│   │   └── {sample}
│   │   │   ├── {sample}_R1.fastq.gz -> 
│   │   │   └── {sample}_R2.fastq.gz -> 
│   └── trimmed
│       ├── {sample}
│       │   ├── {sample}_R1_fastqc.html
│       │   ├── {sample}_R1.fastq.gz
│       │   ├── {sample}_R2_fastqc.html
│       │   └── {sample}_R2.fastq.gz
│       ├── {sample}_R1.fastq.gz_trimming_report.txt
│       ├── {sample}_R1_val_1_fastqc.zip
│       ├── {sample}_R2.fastq.gz_trimming_report.txt
│       └── {sample}_R2_val_2_fastqc.zip
├── logs
├── proactiv
│   ├── plots
│   │   ├── {cell_line}
│   │   ├── promoter_activity_number_hist_all.pdf
│   │   ├── promoter_activity_number_hist_without1.pdf
│   │   ├── promoter_activity_tsne_plot.pdf
│   │   └── tsne_plot_base.pdf
│   └── quantify
│       ├── absolute_promoter_activity.rds
│       ├── alternative_promoters_downReg.rds
│       ├── alternative_promoters_upReg.rds
│       ├── assays_list.rds
│       ├── colData.rds
│       ├── full_alternative_promoters.rds
│       ├── gene_expression.rds
│       ├── goodness_of_fit_ks_pvalues.rds
│       ├── NB_size_per_gene.rds
│       ├── normalized_promoter_counts.rds
│       ├── per_sample_pnbinom_pvalues_matrix.rds
│       ├── proactiv_result.rds
│       ├── promoter_fitted_mu.rds
│       ├── raw_promoter_counts.rds
│       ├── relative_promoter_activity.rds
│       ├── rowData.rds
│       └── sanity_check_DESeq2_result.rds
├── salmon
│   ├── {sample}
│   │   ├── aux_info
│   │   ├── cmd_info.json
│   │   ├── lib_format_counts.json
│   │   ├── libParams
│   │   ├── logs
│   │   └── quant.sf
│   ├── counts
│   │   └── {sample}_promoter_counts.rds
│   ├── merged
│   │   ├── absolute_promoter_activity_category.rds
│   │   ├── absolute_promoter_activity_mean.rds
│   │   ├── absolute_promoter_activity.rds
│   │   ├── dispersion_per_gene.rds
│   │   ├── gene_expression_mean.rds
│   │   ├── gene_expression.rds
│   │   ├── goodness_of_fit_ks_pvalues.rds
│   │   ├── NB_size_per_gene.rds
│   │   ├── normalized_promoter_counts.rds
│   │   ├── per_sample_pnbinom_pvalues_matrix.rds
│   │   ├── promoter_counts.rds
│   │   ├── promoter_fitted_mu.rds
│   │   ├── relative_promoter_activity_mean.rds
│   │   ├── relative_promoter_activity.rds
│   │   ├── sanity_check_DESeq2_result.rds
│   │   └── size_factors.rds
│   └── plots
        ├── {cell_line}
        ├── promoter_activity_number_hist_all.pdf
        ├── promoter_activity_number_hist_without1.pdf
        └── promoter_activity_tsne_plot.pdf
```

- **`benchmarks`** - Benchmarking data for pipeline steps, including runtime and memory usage per sample (e.g., `benchmarks/dexseq_counts.{sample}.txt`).

- **`logs`** - Detailed log files for individual pipeline steps per sample (e.g., `logs/star.{sample}.log`).


# Read accessing and processing
- **`fastqc`** - FastQC quality control reports for raw FASTQ files (e.g., `fastqc/{sample}/{sample}_R1_fastqc.html`,`fastqc/{sample}/{sample}_R1_fastqc.zip`).

- **`fastqs`** - Preparing and processing raw FASTQ files, including:
  - **`raw`** - Original input FASTQ files copied from `--input_dir` (e.g., `raw/{sample}/{sample}_R1.fastq.gz -> /abs/path/to/fastq.gz`).
  - **`downsampled`** - Downsampled FASTQ files after downsampling according to `--downsample_size` (e.g., `trimmed/{sample}/{sample}_R1.fastq.gz`).
  - **`trimmed`** - Processed FASTQ files after initial QC or filtering (e.g., `trimmed/{sample}/{sample}_R1.fastq.gz`). Including FastQC quality control (e.g., `trimmed/{sample}_R1_val_1_fastqc.zip`) and trimming reports (e.g., `trimmed/{sample}_R1.fastq.gz_trimming_report.txt`.)


# Mapping reads
- **`bam`** - Alignment, strand, and juction information, including STAR-aligned BAM file (e.g., `{sample}/{sample}.sorted.bam`),strandness determination for the sample (e.g., `{sample}/{sample}.strand.txt`), and STAR splice junction annotation (e.g., `junctions/{sample}._STARgenome`).


# proActiv
- **`quantify`** - Count promoter reads according to promoter to transcript annotation and only including reads that span first junction, using function of proActiv.
  - Results are normalized and relative activity is generated for comparison (e.g., `normalized_promoter_counts.rds`). 
  - Samples are merged together (e.g., `absolute_promoter_activity.rds`, `assays_list.rds`, `colData.rds`, `proactiv_result.rds`, `rowData.rds`, `gene_expression.rds`). 
  - Sanity check is done by KS test for p value and permutation of condition for false positives. Final numbers are printed in log but complete results may be also accessed (e.g., `promoter_fitted_mu.rds`,`NB_size_per_gene.rds`, `per_sample_pnbinom_pvalues_matrix.rds`, `goodness_of_fit_ks_pvalues.rds`,  `sanity_check_DESeq2_Result.rds`).
  - Alternative promoters are analyzed from absolute and relative promoter activity if the reference condition is specified (e.g., `alternative_promoters_downReg.rds`,`alternative_promoters_downReg.rds`, `full_alternative_promoters.rds`).

## proActiv_plots
- **`plots`** - Visualizes the comparison between promoter activity for different conditions, including 
  - The activity level for different promoter categories (e.g., `{sample}_promoter_activity_category_comparison.pdf`).
  - The number of promoters in each category (e.g., `{sample}_promoter_activity_category_percentage_genewise.pdf`, `sample_promoter_activity_category_percentage.pdf`).
  - The correlation of promoter activity to gene expression for every sample (e.g., `{sample}_promoter_activity_geneexpression_correlation.pdf`).
  - The percentage for each category of promoters at different positions (e.g., `{sample}_activity_position_category.pdf`).
  - The percentage for different categories of gene with promoter(s) (e.g., `{sample}_promoter_activity_single_multiple_category.pdf`).
  - The clustering of samples from the same cell line (e.g., `promoter_activity_tsne_plot.pdf`).
  - The distribution of number of promoters each gene has (e.g., `promoter_activity_number_hist_all.pdf`, `promoter_activity_number_hist_without1.pdf`).

# DEXSeq
- **`{sample}`**
  - Feature count counts the reads falling into each of flattened exon bins (e.g., `{sample}_counts.txt`). 
  - Count counts promoter reads according to exon bin at first exon and promoterId to transcript mapping (e.g., `sample_promoter_counts.rds`). 

## DEXSeq_merge
- **`merge`** - Merge counts from each sample.
  - Samples are merged together (e.g., `promoter_counts.rds`). 
  - Sample counts are normalized to get absolute activity (e.g., `normalized_promoter_counts.rds`, `size_factors.rds`, `promoter_fitted_mu.rds`, `dispersion_per_gene.rds`, `NB_size_per_gene.rds`).
  - Promoter activity and gene expression (e.g., `absolute_promoter_activity.rds`, `gene_expression.rds`, `relative_promoter_activity.rds`).
  - Sanity check is done by KS test for p value and permutation of condition for false positives. Final numbers are printed in log but complete results may be also accessed (e.g.,  `per_sample_pnbinom_pvalues_matrix.rds`, `goodness_of_fit_ks_pvalues.rds`,  `sanity_check_DESeq2_Result.rds`).
  - Categorize promoters for each condition (e.g., `absolute_promoter_activity_mean.rds`, `relative_promoter_activity_mean.rds`, `gene_expression_mean.rds`).
  
- **`plots`** - Plots are the same as [See Plots](#proActiv_plots).

# Salmon
- **`quant`** - Assign read to transcripts and estimate transcript abundancies using k-mer-based quasi-mapping (e.g., `{sample}/quant.sf`).
- **`counts`** - Counts promoter reads according to quasi-mapping and promoterId to transcript mapping (e.g., {sample}_promoter_counts.rds).
- **`merge`** - Result file types are the same as [See Merge](#DEXSeq_merge).
- **`plots`** - Plots are the same as [See Plots](#proActiv_plots).


## Notes
- FASTQ Naming: Files should be named {sample}_R1.fastq.gz and {sample}_R2.fastq.gz.
- Downsampling: Use --downsample_size for uniform sample comparison.
- Pipeline Selection: Choose All to run all tools, or specify a single tool (e.g., Maps) for targeted analysis.
- Custom Parameters: Use --hicpro_params, --hichipper_params, or --hicdc_params to pass tool-specific options.
- Differential Analysis: Provide --samples_comparison for pairwise sample comparisons (if implemented). Provide --reference_condition for reference condition used by proActiv. Provide --method for the method used for differential analysis.

