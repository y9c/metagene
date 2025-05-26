---
layout: default
title: Getting Started
nav_order: 2
---

# Getting Started
{: .no_toc }

This guide will help you get up and running with the Metagene Analysis Package.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Installation

### Requirements

- Python 3.12 or higher
- pip package manager

### Install from PyPI

```bash
pip install metagene
```

### Install from source

```bash
git clone https://github.com/y9c/metagene.git
cd metagene
pip install -e .
```

---

## Basic Usage

### Python API

#### 1. Load Your Data

```python
import metagene

# Load genomic sites from file
sites = metagene.load_sites(
    "sites.tsv", 
    with_header=True,
    meta_col_index=[0, 1, 2]  # chromosome, position, strand
)
```

#### 2. Load Reference Annotations

```python
# Use built-in reference (auto-downloads if needed)
reference = metagene.load_reference("GRCh38")

# Or load custom GTF file
reference = metagene.load_gtf("custom.gtf")
```

#### 3. Run Analysis

```python
# Map sites to transcripts
results = metagene.map_to_transcripts(sites, reference)

# Normalize positions relative to gene structure
gene_bins, gene_stats, gene_splits = metagene.normalize_positions(results, region="all")

# Show summary statistics
print(f"Gene splits - 5'UTR: {gene_splits[0]:.3f}, CDS: {gene_splits[1]:.3f}, 3'UTR: {gene_splits[2]:.3f}")
print(f"Gene statistics - 5'UTR: {gene_stats['5UTR']}, CDS: {gene_stats['CDS']}, 3'UTR: {gene_stats['3UTR']}")
```

#### 4. Generate Plots

```python
# Create metagene plot
metagene.plot_profile(
    gene_bins, 
    gene_splits,
    output_path="metagene_plot.png"
)
```

### Command Line Interface

#### Basic Analysis

```bash
# Run complete analysis
metagene -i sites.tsv -o results.tsv -r GRCh38 -p plot.png --with-header
```

#### Available Options

```bash
# View all options
metagene --help

# List available built-in references
metagene --list

# Download specific reference
metagene --download GRCm39

# Custom analysis parameters
metagene -i sites.bed -o results.tsv -g custom.gtf --bins 200 --region cds
```

---

## Input Data Formats

### Genomic Sites File

Your input file should contain genomic coordinates. Supported formats:

#### BED Format (0-based)
```
chr1    1000    1001    .    100    +
chr1    2000    2001    .    150    -
```

#### TSV Format (1-based with header)
```
chromosome    position    strand    score
chr1          1001        +         100
chr1          2001        -         150
```

#### TSV Format (0-based without header)
```
chr1    1000    1001    +    100
chr1    2000    2001    -    150
```

### Column Specification

Use the `--meta-columns` option to specify which columns contain your coordinates:

```bash
# For BED-like format (chr, start, end, strand)
metagene -i data.tsv -m "1,2,3,4"

# For position-only format (chr, pos, strand)
metagene -i data.tsv -m "1,2,3"
```

---

## Built-in References

The package includes built-in references for common model organisms:

### Available Species

| Species | Reference | Description |
|---------|-----------|-------------|
| **Human** | GRCh38 | Genome Reference Consortium Human Build 38 |
| **Mouse** | GRCm39 | Genome Reference Consortium Mouse Build 39 |
| **Mouse** | GRCm38 | Genome Reference Consortium Mouse Build 38 |
| **Arabidopsis** | TAIR10 | The Arabidopsis Information Resource |
| **Rice** | IRGSP-1.0 | International Rice Genome Sequencing Project |

### Using Built-in References

```python
# List available references
available = metagene.load_reference()
print(available)

# Load specific reference (auto-downloads)
reference = metagene.load_reference("GRCh38")
```

```bash
# Command line
metagene --list
metagene --download GRCh38
```

---

## Examples

### Example 1: RNA Modification Sites

```python
import metagene

# Load m6A modification sites
sites = metagene.load_sites("m6a_sites.bed")

# Use human reference
reference = metagene.load_reference("GRCh38")

# Focus on 3' UTR regions
results = metagene.map_to_transcripts(sites, reference)
gene_bins, gene_stats, gene_splits = metagene.normalize_positions(results, region="3utr")

# Create publication-ready plot
metagene.plot_profile(
    gene_bins, 
    gene_splits,
    "m6a_metagene.png"
)
```

### Example 2: Protein Binding Sites

```bash
# Command line analysis of ChIP-seq peaks
metagene -i chip_peaks.bed \
         -o binding_analysis.tsv \
         -r GRCm39 \
         -p binding_plot.png \
         --region cds \
         --bins 100
```

### Example 3: Custom GTF File

```python
import metagene

# Load sites and custom annotation
sites = metagene.load_sites("sites.tsv", with_header=True)
reference = metagene.load_gtf("custom_annotation.gtf")

# Run analysis
results = metagene.map_to_transcripts(sites, reference)
gene_bins, gene_stats, gene_splits = metagene.normalize_positions(results)

# Generate plot
metagene.plot_profile(gene_bins, gene_splits, "custom_analysis.png")
```

---

## Next Steps

- 📖 **[API Reference](api-reference.md)**: Detailed function documentation
- 🎨 **[Plotting Guide](plotting.md)**: Customizing your visualizations
- 🔧 **[Advanced Usage](advanced.md)**: Power user features
- 💡 **[Examples](examples.md)**: More complete workflows
