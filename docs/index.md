---
layout: default
title: Home
nav_order: 1
description: "Metagene Analysis Package - A Python package for metagene profiling analysis and visualization"
permalink: /
---

# Metagene Analysis Package
{: .fs-9 }

A powerful Python package for metagene profiling analysis and visualization, featuring automatic reference downloading, rich CLI interface, and comprehensive plotting capabilities.
{: .fs-6 .fw-300 }

[Get started now]({{ site.baseurl }}{% link getting-started.md %}){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View it on GitHub](https://github.com/y9c/metagene){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Quick Start

### Installation

```bash
pip install metagene
```

### Basic Usage

```python
import metagene

# Load genomic sites
sites = metagene.load_sites("sites.tsv", with_header=True)

# Load reference annotations
reference = metagene.load_reference("GRCh38")  # Auto-downloads if needed

# Run analysis
results = metagene.map_to_transcripts(sites, reference)
normalized = metagene.normalize_positions(results)

# Generate plot
metagene.simple_metagene_plot(normalized, "output.png")
```

### Command Line Interface

```bash
# Run analysis with built-in reference
metagene -i sites.tsv -o results.tsv -r GRCh38 -p plot.png --with-header

# List available references
metagene --list

# Download specific reference
metagene --download GRCm39
```

---

## Key Features

### 🚀 **Easy Installation & Usage**
One-command installation with pip, intuitive API design, and comprehensive documentation.

### 📊 **Rich Visualization**
Beautiful metagene plots with customizable styling, multiple plot types, and publication-ready output.

### 🧬 **Built-in References**
Automatic downloading of reference genomes for human, mouse, and other model organisms.

### ⚡ **High Performance**
Built on PyRanges v1 and Polars for fast genomic operations and data processing.

### 🔧 **Flexible CLI**
Rich command-line interface with progress bars, emoji feedback, and extensive customization options.

### 🧪 **Well Tested**
Comprehensive test suite with example data and workflows for reliable analysis.

---

## About

Metagene analysis is a powerful technique for visualizing the distribution of genomic features (like protein binding sites, histone modifications, or RNA modifications) relative to gene structures. This package provides a complete toolkit for:

- **Data Loading**: Support for BED, GTF, and TSV formats
- **Reference Management**: Automatic downloading of built-in references
- **Analysis Pipeline**: Complete workflow from raw data to publication-ready plots
- **Customization**: Extensive options for analysis parameters and visualization

### Why Metagene?

- **Simple**: Clean API that's easy to learn and use
- **Fast**: Optimized for large genomic datasets
- **Comprehensive**: Everything you need in one package
- **Reliable**: Well-tested with example workflows
- **Modern**: Built with the latest Python genomics tools

---

## Getting Help

- 📖 **[Documentation]({{ site.baseurl }}{% link getting-started.md %})**: Complete guides and API reference
- 🐛 **[Issues](https://github.com/y9c/metagene/issues)**: Report bugs or request features
- 💬 **[Discussions](https://github.com/y9c/metagene/discussions)**: Ask questions and share ideas
