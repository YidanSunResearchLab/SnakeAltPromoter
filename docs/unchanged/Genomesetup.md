### `docs/Genomesetup.md`

# Genomesetup: Genome Index Preparation

`Genomesetup` sets up genome indices and restriction enzyme fragments for HiChIP/PLAC-seq analysis using Chromap, Hic-pro and MAPS scripts.


## Usage
```bash
Genomesetup [options]
```

### Required Arguments
- `--organism`: Genome assembly (e.g., `hg38`, `mm10`).
- `--organism_fasta`: Absolute path to genome FASTA file (e.g., `/tmp/Genomes/GRCh38/genome.fa`). 

### Optional Arguments
- `-o/--output_dir`: Directory to store genome files (e.g., `~/syidan/Projects/genome`).
- `--restriction_enzyme`: Restriction enzyme (e.g., 'hindiii', 'dpnii', 'bglii', 'ncoi', 'msei', 'hinfI', 'mnase', 'arima'). Default: `mboi`.
- `--bin_size`: Bin size in base pairs (e.g., `5000`). Default: 5000.
- `--threads`: Number of CPU threads (e.g., `30`). Default: 30.

### Example
```bash
Genomesetup \
  --organism hg38 \
  --organism_fasta /path/to/your/genome.fa \
  -o ./genome \
  --restriction_enzyme mboi \
  --bin_size 5000 \
  --threads 30
```

## Output
```
tree organisms/hg38
organisms/hg38
├── BowtieIndex
│   ├── hg38.1.bt2
│   ├── hg38.2.bt2
│   ├── hg38.3.bt2
│   ├── hg38.4.bt2
│   ├── hg38.rev.1.bt2
│   └── hg38.rev.2.bt2
├── ChromapIndex
│   └── hg38.index
├── GenmapIndex
│   ├── index.ids.concat
│   ├── index.ids.limits
│   ├── index.info.concat
│   ├── index.info.limits
│   ├── index.lf.drp
│   ├── index.lf.drp.sbl
│   ├── index.lf.drs
│   ├── index.lf.drv
│   ├── index.lf.drv.sbl
│   ├── index.lf.pst
│   ├── index.rev.lf.drp
│   ├── index.rev.lf.drp.sbl
│   ├── index.rev.lf.drs
│   ├── index.rev.lf.drv
│   ├── index.rev.lf.drv.sbl
│   ├── index.rev.lf.pst
│   ├── index.sa.ind
│   ├── index.sa.len
│   ├── index.sa.val
│   ├── index.txt.concat
│   └── index.txt.limits
├── GenmapIndex.done
├── hg38.autosome.count.txt
├── hg38.chrom.sizes
├── hg38.fa
├── hg38.fa.fai
├── hg38_mappability_50.bedgraph
├── hg38_mappability_50.bw
├── hg38_mappability_50.chrom.sizes
├── hg38_mappability_50.wig
├── main_chroms.txt
└── mboi
    └── 5000
        ├── hg38.txt
        ├── restriction_enzyme.fithichip.bed
        ├── restriction_enzyme.hicdcplus.bintolen.txt.gz
        ├── restriction_enzyme.organism.txt
        ├── restriction_enzyme_position_all.bed
        ├── restriction_enzyme_position_all.map
        ├── restriction_enzyme_position_all.map.filter
        ├── restriction_enzyme_position_all.map.tab
        ├── restriction_enzyme_position_all.txt
        └── restriction_enzyme_resolution_fragments.bed

```
- **`hg38.fa`**: Genome FASTA file.
- **`hg38.chrom.sizes`**: Genome chromosome size file.
- **`BowtieIndex/hg38`**: Bowtie2 index for mapping.
- **`ChromapIndex/hg38.index`**: Chromap index for mapping.
- **`GenmapIndex`**: Genmap index for calculating mappability.
- **`mboi/5000/restriction_enzyme_position_all.txt`**: Restriction fragment coordinates for Maps.
- **`mboi/5000/restriction_enzyme_resolution_fragments.bed`**: Restriction fragment coordinates for Hicpro.
- **`mboi/5000/restriction_enzyme.hicdcplus.bintolen.txt.gz`**: Restriction fragment coordinates for Hicdcplus.
- **`genomesetup_run_hg38_*.log`**: Log file with run details.

## Notes
- Ensure the FASTA file has chromosome names prefixed with `chr` (e.g., `chr1`).
- Use this output folder as `--genome_dir` in `Hichipsnake.py`.
