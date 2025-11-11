# setup.py
from setuptools import setup, find_packages

setup(
    name="SnakeAltPromoter",
    version="1.0.0",
    description="A Snakemake pipeline for alternative promoter analysis",
    author="Yidan Sun",
    author_email="syidan@wustl.edu",
    url="https://github.com/YidanSunResearchLab/SnakeAltPromoter",  # Updated username
    packages=find_packages(),  # Includes workflows/
    python_requires=">=3.10",  # Bumped to 3.10 for modern support
    include_package_data=True,
    install_requires=[
        "snakemake>=8.28.0",  # Core dependency
        "streamlit",          # UI dependencies now default
        "pyarrow",
        "pandas",
    ],

    package_data={
        "snakealtpromoter": [
            "rules/*.Snakefile",       
            "rules/envs/*.yaml",      
            "ui/*.py",
            "scripts/*",
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
        "Programming Language :: Python :: 3.10",  # Matches python_requires
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
