# Metagene

[![Pypi Releases](https://img.shields.io/pypi/v/metagene.svg)](https://pypi.python.org/pypi/metagene)
[![Downloads](https://static.pepy.tech/badge/metagene)](https://pepy.tech/project/metagene)

**Metagene Profiling Analysis and Visualization**

A Python package for performing metagene analysis on genomic sites. This tool allows you to analyze the distribution of genomic features relative to gene regions (5'UTR, CDS, 3'UTR) and create publication-ready metagene profile plots.

## Installation

Install metagene using pip:

```bash
pip install metagene
```

Or using uv:

```bash
uv add metagene
```

## Quick Start

### Command Line Interface

Basic metagene analysis using a built-in reference:

```bash
# Using built-in human genome reference (GRCh38)
metagene -i sites.tsv.gz -r GRCh38 --with-header -m 1,2,3 -w 5 \
         -o output.tsv -s scores.tsv -p plot.png
```

Using a custom GTF file:

```bash
# Using custom GTF annotation
metagene -i sites.bed -g custom.gtf.gz -m 1,2,3 -w 5 \
         -o output.tsv -s scores.tsv -p plot.png
```

### Python API

```python
import polars as pl
from metagene import (
    load_sites, load_reference, map_to_transcripts, 
    normalize_positions, plot_profile
)

# Load your genomic sites
sites_df = load_sites("sites.tsv.gz", with_header=True, meta_col_index=[0, 1, 2])

# Load reference genome annotation
reference = load_reference("GRCh38")  # or load_gtf("custom.gtf.gz")

# Perform metagene analysis
annotated_df = map_to_transcripts(sites_df, reference)
final_df, gene_splits = normalize_positions(annotated_df, strategy="median")

# Generate plot
plot_profile(final_df, gene_splits, "metagene_plot.png")

print(f"Analyzed {len(final_df)} sites")
print(f"Gene splits - 5'UTR: {gene_splits[0]:.3f}, CDS: {gene_splits[1]:.3f}, 3'UTR: {gene_splits[2]:.3f}")
```

## Input Formats

### TSV Format (Tab-separated values)
```
ref	pos	strand	score	pvalue
chr1	1000000	+	0.85	0.001
chr1	2000000	-	0.72	0.005
```

### BED Format
```
chr1	999999	1000000	score1	0.85	+
chr1	1999999	2000000	score2	0.72	-
```

### Column Specification
- Use `-m/--meta-columns` to specify coordinate columns (1-based indexing)
- Use `-w/--weight-columns` to specify score/weight columns
- Use `--with-header` if your file has a header line

## Built-in References

Metagene includes pre-processed gene annotations for major model organisms:

| Species | Assembly | Reference |
|---------|----------|-----------|
| **Human** | GRCh38/hg38 | `GRCh38`, `hg38` |
| | GRCh37/hg19 | `GRCh37`, `hg19` |
| **Mouse** | GRCm39/mm39 | `GRCm39`, `mm39` |
| | GRCm38/mm10 | `GRCm38`, `mm10` |
| | mm9/NCBIM37 | `mm9`, `NCBIM37` |
| **Arabidopsis** | TAIR10 | `TAIR10` |
| **Rice** | IRGSP-1.0 | `IRGSP-1.0` |
| **Model Organisms** | Various | `dm6`, `ce11`, `WBcel235`, `sacCer3`, etc. |

### Managing References

List all available references:
```bash
metagene --list
```

This will show all 23+ available references organized by species:
```
Human:
  GRCh37          - Human genome GRCh37 (Ensembl release 75)
  GRCh38          - Human genome GRCh38 (Ensembl release 110)
  hg19            - Human genome hg19 (UCSC 2021)
  hg38            - Human genome hg38 (UCSC 2022)

Mouse:
  GRCm38          - Mouse genome GRCm38 (Ensembl release 102)
  GRCm39          - Mouse genome GRCm39 (Ensembl release 110)
  mm10            - Mouse genome mm10 (UCSC 2021)
  mm39            - Mouse genome mm39 (UCSC 2024)
  mm9             - Mouse genome mm9 (UCSC 2020)

... and more
```

Download a specific reference:
```bash
metagene --download GRCh38
```

Download all references (requires ~10GB disk space):
```bash
metagene --download all
```


## CLI Examples

### Basic Analysis

```bash
# Analyze sites with built-in human reference
metagene -i sites.tsv.gz -r GRCh38 --with-header \
         -m 1,2,3 -w 5 -o output.tsv -p plot.png
```

Note: References are automatically downloaded on first use.

### Advanced Options
```bash
# Full analysis with custom parameters
metagene -i sites.bed -r GRCh38 \
         -m 1,2,3 -w 5,6 -n "score1,score2" \
         --bins 200 --region all \
         --score-transform log2 --normalize \
         -o annotated.tsv -s statistics.tsv -p metagene.pdf
```

### Custom GTF Reference
```bash
# Use your own GTF annotation
metagene -i sites.tsv.gz -g annotation.gtf.gz --with-header \
         -m 1,2,3 -w 4 -o output.tsv -p plot.png
```

## API Reference

### Core Functions

- `load_sites(file, with_header=False, meta_col_index=[0,1,2])` - Load genomic sites
- `load_reference(name)` - Load built-in reference genome
- `load_gtf(file)` - Load custom GTF annotation  
- `map_to_transcripts(sites, reference)` - Annotate sites with gene information
- `normalize_positions(annotated_sites, strategy="median")` - Normalize to relative positions
- `plot_profile(data, gene_splits, output_file)` - Generate metagene plot

### Analysis Workflow

```python
# 1. Load data
sites = load_sites("input.tsv", with_header=True, meta_col_index=[0,1,2])
reference = load_reference("GRCh38")

# 2. Annotate and normalize  
annotated = map_to_transcripts(sites, reference)
normalized, splits = normalize_positions(annotated)

# 3. Visualize
plot_profile(normalized, splits, "output.png")
```

## Demo

![Metagene Profile](docs/fig_metagene.svg)

The plot shows the distribution of genomic sites across normalized gene regions:
- **5'UTR** (0.0 - first split): 5' untranslated region
- **CDS** (first split - second split): Coding sequence  
- **3'UTR** (second split - 1.0): 3' untranslated region