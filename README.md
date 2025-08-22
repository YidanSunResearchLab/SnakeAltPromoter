# SnakeAltPromoter Facilitates Differential Alternative Promoter Analysis

**SnakeAltPromoter** is a Snakemake-based pipeline designed for end-to-end analysis of alternative promoter. It streamlines the processing of RNA-seq data by integrating **FastQC**, **TrimGlore**, **STAR** and **multiQC**, and supports multiple state-of-the-art promoter counting algorithms including **Proactiv**, **DEXseq** and **Salmon**. Finally it includes differential analysis using **DESeq2** and **Proactiv**. This pipeline offers researchers a reproducible, modular, and scalable solution for analyzing alternative promoter for RNA-seq data across diverse experimental designs.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/YidanSunResearchLab/SnakeAltPromoter.git
   cd SnakeAltPromoter
   ```

2. **Install Dependencies**:
   - Set up the environment:
  ```bash
  conda create -n SnakeAltPromoter -c conda-forge python>=3.8
  conda activate SnakeAltPromoter
  pip install .
  ```

3. **Verify Installation**:
   ```bash
   Snakealtpromoter --help
   ```

## Quick Start

SnakeAltPromoter offers three main commands:

### 1. Genome Setup
Prepare genome indices and promoter annotation:
```bash
Genomesetup \
  --organism hg38 \
  --organism_fasta /path/to/your/genome.fa \
  --genes_gtf /path/to/your/genes.gtf \
  -o ./genome \
  --threads 30
```

### 2.  RNA-Seq Processing for alternative promoter analysis
Process FASTQs:
```bash
Snakealtpromoter.py \
  -i /path/to/input/fastqs/dir/ 
  --genome_dir /path/to/genomesetup/dir/ 
  -o /path/to/output/dir/ 
  --threads 30 
  --organism hg38 
  --control_condition Healthy 
  --test_condition Heart_Failure 
```

## Documentation
For detailed usage, see:
- [Genomesetup.py](docs/Genomesetup.md)
- [Snakealtpromoter.py](docs/Snakealtpromoter.md)

## Contributing
Feel free to open issues or submit pull requests on [GitHub](https://github.com/YidanSunResearchLab/SnakeAltPromoter).

## Citation
Please cite this paper: SnakeAltPromoter Facilitates Differential Alternative Promoter Analysis. [BioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.16.669128v1).

## License
See [LICENSE](LICENSE) for details.


