---
layout: default
title: Home
nav_order: 1
description: 'Metagene Analysis Package - A Python package for metagene profiling analysis and visualization'
permalink: /
---

# Metagene Analysis Package

{: .fs-9 }

A powerful Python package for metagene profiling analysis and visualization, featuring automatic reference downloading, rich CLI interface, and comprehensive plotting capabilities.
{: .fs-6 .fw-300 }

[Get started now](getting-started.md){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
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
metagene.plot_profile(normalized, "output.png")
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

## Getting Help

- 📖 **[Documentation](getting-started.md)**: Complete guides and API reference
- 🐛 **[Issues](https://github.com/y9c/metagene/issues)**: Report bugs or request features
- 💬 **[Discussions](https://github.com/y9c/metagene/discussions)**: Ask questions and share ideas
