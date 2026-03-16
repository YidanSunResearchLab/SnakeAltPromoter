# Validate sample sheet error handling
import os
import subprocess
import sys

from test_helpers import (
    create_fastqs,
    create_fake_snakemake,
    create_genome_dir,
    env_with_fake_snakemake,
)


def run_snakealtpromoter(args, env):
    cmd = [
        sys.executable,
        "-m",
        "snakealtpromoter.workflows.Snakealtpromoter"
    ] + args
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


def test_samplesheet_missing_columns(tmp_path):
    bin_dir = tmp_path / "bin"
    create_fake_snakemake(bin_dir)
    env = env_with_fake_snakemake(os.environ, bin_dir)

    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    genome_dir = create_genome_dir(tmp_path)
    create_fastqs(input_dir, paired=True)

    bad_sheet = tmp_path / "bad_samplesheet.tsv"
    bad_sheet.write_text("sampleName\tcondition\nsample1\toverall\n")

    args = [
        "-i",
        str(input_dir),
        "-o",
        str(output_dir),
        "--organism",
        "hg38",
        "--genome_dir",
        str(genome_dir),
        "--reads",
        "paired",
        "--method",
        "salmon",
        "--sample_sheet",
        str(bad_sheet),
        "--threads",
        "2",
    ]

    result = run_snakealtpromoter(args, env)
    assert result.returncode != 0
    assert "missing columns" in result.stderr


def test_samplesheet_missing_sample(tmp_path):
    bin_dir = tmp_path / "bin"
    create_fake_snakemake(bin_dir)
    env = env_with_fake_snakemake(os.environ, bin_dir)

    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    genome_dir = create_genome_dir(tmp_path)
    create_fastqs(input_dir, paired=True)

    sheet = tmp_path / "samplesheet.tsv"
    sheet.write_text(
        "sampleName\tcondition\tbatch\tdifferential\n"
        "not_in_inputs\toverall\tbatch1\tcontrol\n"
    )

    args = [
        "-i",
        str(input_dir),
        "-o",
        str(output_dir),
        "--organism",
        "hg38",
        "--genome_dir",
        str(genome_dir),
        "--reads",
        "paired",
        "--method",
        "salmon",
        "--sample_sheet",
        str(sheet),
        "--threads",
        "2",
    ]

    result = run_snakealtpromoter(args, env)
    assert result.returncode != 0
    assert "Samples in sampleSheet but not in input_dir" in result.stderr


def test_samplesheet_valid_two_conditions(tmp_path):
    bin_dir = tmp_path / "bin"
    create_fake_snakemake(bin_dir)
    env = env_with_fake_snakemake(os.environ, bin_dir)

    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    genome_dir = create_genome_dir(tmp_path)
    create_fastqs(input_dir, sample="sample1", paired=True)
    create_fastqs(input_dir, sample="sample2", paired=True)

    sheet = tmp_path / "samplesheet.tsv"
    sheet.write_text(
        "sampleName\tcondition\tbatch\tdifferential\n"
        "sample1\tcontrol_cond\tbatch1\tcontrol\n"
        "sample2\ttest_cond\tbatch1\ttest\n"
    )

    args = [
        "-i",
        str(input_dir),
        "-o",
        str(output_dir),
        "--organism",
        "hg38",
        "--genome_dir",
        str(genome_dir),
        "--reads",
        "paired",
        "--method",
        "salmon",
        "--sample_sheet",
        str(sheet),
        "--threads",
        "2",
    ]

    result = run_snakealtpromoter(args, env)
    assert result.returncode == 0, result.stderr
    assert (output_dir / "Snakealtpromoter_config.yaml").exists()
