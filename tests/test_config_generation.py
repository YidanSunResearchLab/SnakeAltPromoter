# Validate config creation for all methods and reads combinations
# without running the real snakemake
import os
import subprocess
import sys

import yaml

from test_helpers import (
    create_fastqs,
    create_fake_snakemake,
    create_genome_dir,
    env_with_fake_snakemake,
)


def write_samplesheet(
    path, sample="sample1", condition="overall", batch="batch1", role="control"
):
    path.write_text(
        "\t".join(["sampleName", "condition", "batch", "differential"])
        + "\n"
        + "\t".join([sample, condition, batch, role])
        + "\n"
    )


def run_snakealtpromoter(args, env):
    cmd = [
        sys.executable,
        "-m",
        "snakealtpromoter.workflows.Snakealtpromoter"
    ] + args
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


def test_config_generation_all_methods(tmp_path):
    bin_dir = tmp_path / "bin"
    create_fake_snakemake(bin_dir)
    env = env_with_fake_snakemake(os.environ, bin_dir)

    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    sample_sheet = tmp_path / "samplesheet.tsv"
    genome_dir = create_genome_dir(tmp_path)

    for method in ["salmon", "proactiv", "dexseq", "cage", "rnaseq"]:
        for reads in ["single", "paired"]:
            input_dir_method = input_dir / f"{method}_{reads}"
            output_dir_method = output_dir / f"{method}_{reads}"
            input_dir_method.mkdir(parents=True, exist_ok=True)
            output_dir_method.mkdir(parents=True, exist_ok=True)

            create_fastqs(input_dir_method, paired=(reads == "paired"))
            write_samplesheet(sample_sheet)

            args = [
                "-i",
                str(input_dir_method),
                "-o",
                str(output_dir_method),
                "--organism",
                "hg38",
                "--genome_dir",
                str(genome_dir),
                "--reads",
                reads,
                "--method",
                method,
                "--sample_sheet",
                str(sample_sheet),
                "--threads",
                "2",
            ]

            result = run_snakealtpromoter(args, env)
            assert result.returncode == 0, result.stderr

            config_path = output_dir_method / "Snakealtpromoter_config.yaml"
            assert config_path.exists()
            config = yaml.safe_load(config_path.read_text())

            assert config["method"] == method
            if reads == "paired":
                assert config["reads"] == ["R1", "R2"]
            else:
                assert config["reads"] == ["R1"]
