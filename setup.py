
# setup.py
from setuptools import setup, find_packages

setup(
    name="SnakeAltPromoter",
    version="0.1.0",
    description="A Snakemake pipeline for alternative promoter analysis",
    author="Yidan Sun",
    author_email="syidan@wustl.edu",
    url="https://github.com/YidanSunResearchLab/SnakeAltPromoter",  # Updated username
    packages=find_packages(),  # Includes workflows/
    install_requires=[
        "snakemake>=8.28.0",  # Core dependency
    ],
    python_requires=">=3.8",  # Bumped to 3.8 for modern support
    package_data={
        "workflows": [
            "../scripts/*", 
            "../rules/*.Snakefile",       # Snakemake workflows
            "../rules/envs/*.yaml",       # Conda envs
            #"../organisms/*"              # Small organism files only
        ],
    },
    entry_points={
        "console_scripts": [
            "Genomesetup = workflows.Genomesetup:main",
            "Snakealtpromoter = workflows.Snakealtpromoter:main",
        ]
    },
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",  # Matches python_requires
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
