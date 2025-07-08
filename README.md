# OrbitAltPromoter: A Comprehensive Pipeline for Alternative Promoter Analysis

**OrbitAltPromoter** is a Snakemake-based pipeline designed for end-to-end analysis of alternative promoter. It streamlines the processing of RNA-seq data by integrating robust mapping tools such as **STAR** and **Salmon**, and supports multiple state-of-the-art promoter counting algorithms including **Proactiv**, **DEXseq** and **Salmon**. 

OrbitAltPromoter also includes features for **genome setup**, **read downsampling**, and **comprehensive quality control**, offering researchers a reproducible, modular, and scalable solution for analyzing data across diverse experimental designs.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/YidanSunResearchLab/OrbitAltPromoter.git
   cd OrbitAltPromoter
   ```

2. **Install Dependencies**:
   - Set up the environment:
  ```bash
  conda create -n OrbitAltPromoter -c conda-forge python>=3.8
  conda activate OrbitAltPromoter
  pip install .
  ```

3. **Verify Installation**:
   ```bash
   Altpromoterflow --help
   ```

## Quick Start

OrbitAltPromoter offers three main commands:

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
Altpromoterflow \
  -i /home/syidan/syidan/Projects/hichip_real_samples/ \
  --genome_dir /abs/path/to/genome \
  -o /home/syidan/syidan/Projects/hichip_real_samples_output/ \
  --threads 100 \
  --organism hg38 \
  --bin_size 5000 \
  --downsample_size 50000000
```

## Documentation
For detailed usage, see:
- [Genomesetup.py](docs/Genomesetup.md)
- [Altpromoterflow.py](docs/Altpromoterflow.md)

## Contributing
Feel free to open issues or submit pull requests on [GitHub](https://github.com/YidanSunResearchLab/HiChIPSnake).

## License
See [LICENSE](LICENSE) for details.


