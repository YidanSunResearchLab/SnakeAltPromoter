from importlib import resources


def test_data_files_accessible():
    data_root = resources.files("snakealtpromoter") / "data"

    expected = [
        data_root / "hg38.chr.txt",
        data_root / "hg38.fa",
        data_root / "hg38.gtf",
        data_root / "samplesheet" / "samplesheet.tsv",
    ]

    for path in expected:
        assert path.is_file(), f"Missing data file: {path}"


def test_docs_files_accessible():
    docs_root = resources.files("snakealtpromoter") / "docs"

    expected = [
        docs_root / "Genomesetup.md",
        docs_root / "Snakealtpromoter.md",
        docs_root / "workflow_overview.png",
    ]

    for path in expected:
        assert path.is_file(), f"Missing docs file: {path}"
