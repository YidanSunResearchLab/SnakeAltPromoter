# Some common utilities that are used by the other files for
# mocking `snakemake` existence and validating processing
# of workflow configs.
import os
import stat


def create_fake_snakemake(bin_dir):
    bin_dir.mkdir(parents=True, exist_ok=True)
    if os.name == "nt":
        script_path = bin_dir / "snakemake.cmd"
        script_path.write_text("@echo off\r\necho fake snakemake\r\nexit /b 0\r\n")
        return script_path

    script_path = bin_dir / "snakemake"
    script_path.write_text("#!/bin/sh\necho 'fake snakemake'\nexit 0\n")
    script_path.chmod(script_path.stat().st_mode | stat.S_IEXEC)
    return script_path


def create_genome_dir(base_dir, organism="hg38"):
    genome_dir = base_dir / "genome"
    ann_dir = genome_dir / "organisms" / organism / "Annotation"
    star_dir = genome_dir / "organisms" / organism / "STARIndex"
    ann_dir.mkdir(parents=True, exist_ok=True)
    star_dir.mkdir(parents=True, exist_ok=True)

    required_files = [
        genome_dir / "organisms" / organism / "Annotation.gtf",
        ann_dir / "genes.bed",
        ann_dir / "DEXSeq_flattened_exons.gff",
        ann_dir / "proActiv_promoter_annotation.rds",
        ann_dir / "genes_t2g.tsv",
    ]
    for path in required_files:
        path.write_text("")

    return genome_dir


def create_fastqs(input_dir, sample="sample1", paired=True):
    input_dir.mkdir(parents=True, exist_ok=True)
    (input_dir / f"{sample}_R1.fastq.gz").write_text("")
    if paired:
        (input_dir / f"{sample}_R2.fastq.gz").write_text("")


def env_with_fake_snakemake(base_env, bin_dir):
    env = dict(base_env)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    return env
