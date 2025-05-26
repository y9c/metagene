#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Configuration data for Metagene

from typing import TypedDict

class ReferenceInfo(TypedDict):
    parquet_file: str
    source_file: str
    description: str

# GitHub repository information
GITHUB_REPO = "y9c/metagene"
GITHUB_API_URL = f"https://api.github.com/repos/{GITHUB_REPO}/releases/latest"
GITHUB_DOWNLOAD_BASE = f"https://github.com/{GITHUB_REPO}/releases/download/data"

# Built-in reference mappings with source information
BUILTIN_REFERENCES: dict[str, ReferenceInfo] = {
    # Human genomes
    "GRCh38": {
        "parquet_file": "GRCh38.parquet",
        "source_file": "Homo_sapiens/raw/GRCh38.release110.gtf.gz",
        "description": "Human genome GRCh38 (Ensembl release 110)"
    },
    "hg38": {
        "parquet_file": "hg38.parquet",
        "source_file": "Homo_sapiens/raw/hg38.20221028.gtf.gz",
        "description": "Human genome hg38 (UCSC 2022)"
    },
    "GRCh37": {
        "parquet_file": "GRCh37.parquet",
        "source_file": "Homo_sapiens/raw/GRCh37.release75.gtf.gz",
        "description": "Human genome GRCh37 (Ensembl release 75)"
    },
    "hg19": {
        "parquet_file": "hg19.parquet",
        "source_file": "Homo_sapiens/raw/hg19.20210517.gtf.gz",
        "description": "Human genome hg19 (UCSC 2021)"
    },
    # Mouse genomes
    "GRCm39": {
        "parquet_file": "GRCm39.parquet",
        "source_file": "Mus_musculus/raw/GRCm39.release110.gtf.gz",
        "description": "Mouse genome GRCm39 (Ensembl release 110)"
    },
    "mm39": {
        "parquet_file": "mm39.parquet",
        "source_file": "Mus_musculus/raw/mm39.20240214.gtf.gz",
        "description": "Mouse genome mm39 (UCSC 2024)"
    },
    "GRCm38": {
        "parquet_file": "GRCm38.parquet",
        "source_file": "Mus_musculus/raw/GRCm38.release102.gtf.gz",
        "description": "Mouse genome GRCm38 (Ensembl release 102)"
    },
    "mm10": {
        "parquet_file": "mm10.parquet",
        "source_file": "Mus_musculus/raw/mm10.20210423.gtf.gz",
        "description": "Mouse genome mm10 (UCSC 2021)"
    },
    "mm9": {
        "parquet_file": "mm9.parquet",
        "source_file": "Mus_musculus/raw/mm9.20200110.gtf.gz",
        "description": "Mouse genome mm9 (UCSC 2020)"
    },
    "NCBIM37": {
        "parquet_file": "NCBIM37.parquet",
        "source_file": "Mus_musculus/raw/NCBIM37.release67.gtf.gz",
        "description": "Mouse genome NCBIM37 (Ensembl release 67)"
    },
    # Model organisms
    "TAIR10": {
        "parquet_file": "TAIR10.parquet",
        "source_file": "Arabidopsis_thaliana/raw/TAIR10.release57.gtf.gz",
        "description": "Arabidopsis thaliana TAIR10"
    },
    "IRGSP-1.0": {
        "parquet_file": "IRGSP-1.0.parquet",
        "source_file": "Oryza_sativa/raw/IRGSP-1.0.release57.gtf.gz",
        "description": "Rice IRGSP-1.0"
    },
    "WBcel235": {
        "parquet_file": "WBcel235.parquet",
        "source_file": "Caenorhabditis_elegans/raw/WBcel235.release110.gtf.gz",
        "description": "C. elegans WBcel235"
    },
    "ce11": {
        "parquet_file": "ce11.parquet",
        "source_file": "Caenorhabditis_elegans/raw/ce11.20200110.gtf.gz",
        "description": "C. elegans ce11 (UCSC)"
    },
    "BDGP6.32": {
        "parquet_file": "BDGP6.32.parquet",
        "source_file": "Drosophila_melanogaster/raw/BDGP6.32.release110.gtf.gz",
        "description": "D. melanogaster BDGP6.32"
    },
    "dm6": {
        "parquet_file": "dm6.parquet",
        "source_file": "Drosophila_melanogaster/raw/dm6.20210210.gtf.gz",
        "description": "D. melanogaster dm6 (UCSC)"
    },
    "GRCz11": {
        "parquet_file": "GRCz11.parquet",
        "source_file": "Danio_rerio/raw/GRCz11.release110.gtf.gz",
        "description": "Zebrafish GRCz11"
    },
    "GRCz10": {
        "parquet_file": "GRCz10.parquet",
        "source_file": "Danio_rerio/raw/GRCz10.release91.gtf.gz",
        "description": "Zebrafish GRCz10"
    },
    "bGalGal1": {
        "parquet_file": "bGalGal1.parquet",
        "source_file": "Gallus_gallus/raw/bGalGal1.mat.broiler.GRCg7b.release110.gtf.gz",
        "description": "Chicken bGalGal1.mat.broiler.GRCg7b"
    },
    "Glycine_max_v2.1": {
        "parquet_file": "Glycine_max_v2.1.parquet",
        "source_file": "Glycine_max/raw/Glycine_max_v2.1.release57.gtf.gz",
        "description": "Soybean Glycine max v2.1"
    },
    "R64-1-1": {
        "parquet_file": "R64-1-1.parquet",
        "source_file": "Saccharomyces_cerevisiae/raw/R64-1-1.release57.gtf.gz",
        "description": "S. cerevisiae R64-1-1"
    },
    "sacCer3": {
        "parquet_file": "sacCer3.parquet",
        "source_file": "Saccharomyces_cerevisiae/raw/sacCer3.20210210.gtf.gz",
        "description": "S. cerevisiae sacCer3 (UCSC)"
    },
    "ASM294v2": {
        "parquet_file": "ASM294v2.parquet",
        "source_file": "Schizosaccharomyces_pombe/raw/ASM294v2.release57.gtf.gz",
        "description": "S. pombe ASM294v2"
    }
} 