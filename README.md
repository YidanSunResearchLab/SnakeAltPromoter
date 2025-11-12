# 🧬 SnakeAltPromoter Facilitates Differential Alternative Promoter Analysis

**SnakeAltPromoter** is a **Snakemake-based pipeline** for streamlined, reproducible, and scalable analysis of **alternative promoter usage** from RNA-seq or CAGE data.
It integrates all major steps—from raw read preprocessing to promoter-level quantification and differential analysis—using state-of-the-art tools.

## 🧭 Workflow Overview

Below is the schematic overview of the **SnakeAltPromoter** pipeline, as shown in our paper:

![Workflow Overview](https://raw.githubusercontent.com/YidanSunResearchLab/SnakeAltPromoter/main/snakealtpromoter/docs/workflow_overview.png)

*Figure 1. Overview of the SnakeAltPromoter pipeline, showing genome setup, RNA-seq/CAGE processing, promoter quantification, classification, and differential promoter usage analysis.*

---

## ✳️ Key Features

* **End-to-end automation** of alternative promoter analysis from FASTQ to results.
* **Integrated QC and preprocessing** via **FastQC**, **TrimGalore**, **STAR**, and **MultiQC**.
* **Supports multiple promoter quantification tools:**

  * [ProActiv](https://github.com/GoekeLab/proActiv)
  * [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html)
  * [Salmon](https://salmon.readthedocs.io/en/latest/)
* **Built-in differential analysis** using **DESeq2** and ProActiv modules.
* **Fully modular, reproducible, and scalable**—ideal for large multi-sample RNA-seq datasets.

---

## 🧩 Installation

### 1. Install via Conda (Recommended)

```bash
# Create and activate environment
conda create -n SnakeAltPromoter -c bioconda SnakeAltPromoter
conda activate SnakeAltPromoter
```

### 2. Install via Pip

```bash
pip install --user --no-cache-dir SnakeAltPromoter
```

### 3. Manual Installation (Build from Source)

```bash
git clone https://github.com/YidanSunResearchLab/SnakeAltPromoter.git
cd SnakeAltPromoter
conda create -n SnakeAltPromoter -c conda-forge python
conda activate SnakeAltPromoter
pip install .
```

### 4. Verify Installation

```bash
Snakealtpromoter --help
```

If successful, usage instructions and command-line options will be displayed.

---

## 🖥️ Optional: Launch GUI

For a graphical interface:

```bash
sap-ui
```

---

## 🧪 Minimal Test Case

Example human chr22 genome fasta and gtf files are available in the `snakealtpromoter/data/` directory.
Example RNA-seq data are available in the `snakealtpromoter/data/` directory.

### 1. Genome Setup

```bash
git clone https://github.com/YidanSunResearchLab/SnakeAltPromoter.git
cd SnakeAltPromoter
Genomesetup \
  --organism hg38 \
  --organism_fasta "$(pwd)/snakealtpromoter/data/hg38.fa" \
  --genes_gtf "$(pwd)/snakealtpromoter/data/hg38.gtf" \
  -o ./genome \
  --threads 30
```

### 2. Run Test Analysis

```bash
git clone https://github.com/YidanSunResearchLab/SnakeAltPromoter.git
cd SnakeAltPromoter
Snakealtpromoter \
  -i "$(pwd)/snakealtpromoter/data" \
  --genome_dir "$(pwd)/genome" \
  -o test_output \
  --threads 30 \
  --organism hg38 --trim \
  --sample_sheet "$(pwd)/snakealtpromoter/data/samplesheet/samplesheet.tsv"
```

Output directory structure is described in the documentation:

* [Genomesetup](snakealtpromoter/docs/Genomesetup.md)
* [Snakealtpromoter](snakealtpromoter/docs/Snakealtpromoter.md)

---

## 📊 Reproduce Results from the Paper

To reproduce analyses in the **SnakeAltPromoter** manuscript:

1. Download the **GENCODE v46** genome FASTA and GTF from [GENCODE](https://www.gencodegenes.org/human/release_46.html).
2. Retrieve **heart RNA-seq** and **CAGE** data from
   [GSE147236](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147236) using `fastq-dump`.

### 1. Genome Setup

```bash
Genomesetup \
  --organism hg38 \
  --organism_fasta /Absolute/path/to/genome.fa \
  --genes_gtf /Absolute/path/to/genes.gtf \
  -o /Absolute/path/to/genome \
  --threads 30
```

### 2. RNA-seq Processing

```bash
git clone https://github.com/YidanSunResearchLab/SnakeAltPromoter.git
cd SnakeAltPromoter
Snakealtpromoter \
  -i /Absolute/path/to/heart_RNAseq/ \
  --genome_dir /Absolute/path/to/genome \
  -o heart_RNAseq_output \
  --threads 30 \
  --organism hg38 --trim \
  --sample_sheet "$(pwd)/snakealtpromoter/data/samplesheet/Heart_RNAseq.tsv"
```

### 3. CAGE Processing

```bash
git clone https://github.com/YidanSunResearchLab/SnakeAltPromoter.git
cd SnakeAltPromoter
Snakealtpromoter \
  -i /Absolute/path/to/heart_CAGE/ \
  --genome_dir /Absolute/path/to/genome \
  -o heart_CAGE_output \
  --threads 30 \
  --organism hg38 \
  --sample_sheet "$(pwd)/snakealtpromoter/data/samplesheet/Heart_CAGE.tsv" \
  --method cage --reads single
```

### 4. Compare with Published Results

| Supplementary Table | Description                                                                  |
| ------------------- | ---------------------------------------------------------------------------- |
| **Table 1**         | Comprehensive promoter coordinates                                           |
| **Table 2**         | Promoter classifications (major/minor) by ProActiv, Salmon, DEXSeq, and CAGE |
| **Table 3**         | Promoter counts across samples                                               |
| **Table 4**         | Differential promoter activity (healthy vs. failure) across tools            |

---

## 🚀 Quick Start for your own genome and sequencing data

### Step 1. Genome Setup

Prepare genome indices and promoter annotations:

```bash
Genomesetup \
  --organism hg38 \
  --organism_fasta /Absolute/path/to/genome.fa \
  --genes_gtf /Absolute/path/to/genes.gtf \
  -o /Absolute/path/to/genome \
  --threads 30
```

### Step 2. Process RNA-seq Data

Run alternative promoter analysis:

```bash
Snakealtpromoter \
  -i /Absolute/path/to/input_fastqs/ \
  --genome_dir /Absolute/path/to/genome \
  -o ./output/ \
  --threads 30 \
  --organism hg38 --trim \
  --sample_sheet /Absolute/path/to/samplesheet.tsv \
  --method cage --reads single   # Add these only for CAGE data
```

For detailed documentation, see:

* [Genomesetup](snakealtpromoter/docs/Genomesetup.md)
* [Snakealtpromoter](snakealtpromoter/docs/Snakealtpromoter.md)

---

## 🤝 Contributing

Contributions are welcome!
Please open an issue or submit a pull request via [GitHub](https://github.com/YidanSunResearchLab/SnakeAltPromoter).

---

## 🧾 Citation

If you use **SnakeAltPromoter**, please cite:

> Tan J. *et al.* (2025). **SnakeAltPromoter Facilitates Differential Alternative Promoter Analysis.**
> *bioRxiv*. [https://doi.org/10.1101/2025.08.16.669128](https://www.biorxiv.org/content/10.1101/2025.08.16.669128v1)

---

## ⚖️ License

See the [LICENSE](LICENSE.md) file for details.
