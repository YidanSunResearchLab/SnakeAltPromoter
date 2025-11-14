# setup.py
from setuptools import setup, find_packages

setup(
    name="SnakeAltPromoter",
    version="1.0.5",
    description="A Snakemake pipeline for alternative promoter analysis",
    author="Yidan Sun",
    author_email="syidan@wustl.edu",
    url="https://github.com/YidanSunResearchLab/SnakeAltPromoter",  # Updated username
    packages=find_packages(),  # Includes workflows/
    python_requires=">=3.11",  # Bumped to 3.11 for modern support
    include_package_data=True,
    install_requires=[
        "snakemake>=8.28.0",  # Core dependency
        "streamlit",          # UI dependencies now default
        "pyarrow",
        "pandas",
    ],

    # When `pyproject.toml` is present, `setup.py` is not the entry point
    # during install. However, it may be still read indirectly for consultation
    #  - i.e. if it's called by `pyproject.toml`. For e.g. all build actions in
    # the README.md - e.g. `pip install .` will use `pyproject.toml` - `setup.py`
    # enables legacy builds via `python setup.py install`. (See PEP 517 and 518)
    package_data={
        "snakealtpromoter": [
            "data/**",
            "docs/**",
            "rules/*.Snakefile",
            "rules/envs/*.yaml",
            "scripts/downstream_plots/**",
        ],
    },
    entry_points={
        "console_scripts": [
            "Genomesetup = snakealtpromoter.workflows.Genomesetup:main",
            "Snakealtpromoter = snakealtpromoter.workflows.Snakealtpromoter:main",
            "sap=snakealtpromoter.cli:main",
            "sap-ui=snakealtpromoter.ui.launch:main",
        ]
    },
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.11",  # Matches python_requires
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
