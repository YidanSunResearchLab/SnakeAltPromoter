# SnakeAltPromoter: A Comprehensive Pipeline for Differential Alternative Promoter Analysis

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
Prepare genome indices and restriction fragments:
```bash
Genomesetup \
  --organism hg38 \
  --organism_fasta /path/to/your/genome.fa \
  -o ./genome \
  --threads 100
```

### 2.  RNA-Seq Processing for alternative promoter analysis
Process FASTQs into interaction maps:
```bash
Snakealtpromoter \
  -i /home/Projects/hichip_real_samples/ \
  --genome_dir /abs/path/to/genome \
  -o /home/Projects/hichip_real_samples_output/ \
  --threads 100 \
  --organism hg38 \
  --bin_size 5000 \
  --downsample_size 50000000
```

## Documentation
For detailed usage, see:
- [Genomesetup.py](docs/Genomesetup.md)
- [Snakealtpromoter.py](docs/Snakealtpromoter.md)

## Contributing
Feel free to open issues or submit pull requests on [GitHub](https://github.com/YidanSunResearchLab/SnakeAltPromoter).

## License
See [LICENSE](LICENSE) for details.


