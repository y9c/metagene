---
layout: default
title: API Reference
nav_order: 3
---

# API Reference
{: .no_toc }

Complete reference for all functions and classes in the Metagene package.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Data I/O Functions

### load_sites()

Load genomic sites from various file formats.

```python
metagene.load_sites(
    input_file_name: str,
    with_header: bool = False,
    meta_col_index: List[int] = None
) -> PyRanges
```

**Parameters:**
- `input_file_name`: Path to input file (BED, TSV, or compressed formats)
- `with_header`: Whether the file has a header row
- `meta_col_index`: Column indices for [chr, start, end, strand] or [chr, pos, strand]

**Returns:**
- PyRanges object with genomic intervals

**Example:**
```python
# Load BED file without header
sites = metagene.load_sites("sites.bed", meta_col_index=[0, 1, 2, 5])

# Load TSV with header
sites = metagene.load_sites("sites.tsv", with_header=True, meta_col_index=[0, 1, 3])
```

### load_reference()

Load built-in reference annotations or list available references.

```python
metagene.load_reference(species: Optional[str] = None) -> Union[PyRanges, dict]
```

**Parameters:**
- `species`: Species name (e.g., "GRCh38", "GRCm39") or None to list available

**Returns:**
- PyRanges with annotations if species specified, dict of available species otherwise

**Example:**
```python
# List available references
available = metagene.load_reference()

# Load specific reference
reference = metagene.load_reference("GRCh38")
```

### load_gtf()

Load custom GTF/GFF file.

```python
metagene.load_gtf(gtf_file: str) -> PyRanges
```

**Parameters:**
- `gtf_file`: Path to GTF or GFF file

**Returns:**
- PyRanges object with gene annotations

---

## Analysis Functions

### map_to_transcripts()

Map genomic sites to transcript coordinates.

```python
metagene.map_to_transcripts(
    input_sites: pr.PyRanges,
    exon_ref: pr.PyRanges
) -> pl.DataFrame
```

**Parameters:**
- `input_sites`: PyRanges with genomic sites (output from load_sites())
- `exon_ref`: PyRanges with gene/transcript annotations (output from load_reference() or load_gtf())

**Returns:**
- Polars DataFrame with transcript mapping information including transcript coordinates

**Example:**
```python
# Load data
sites = metagene.load_sites("sites.tsv", with_header=True, meta_col_index=[0, 1, 2])
reference = metagene.load_reference("GRCh38")

# Map to transcripts  
annotated = metagene.map_to_transcripts(sites, reference)
```

### normalize_positions()

Normalize transcript positions to relative feature positions (0-1 scale).

```python
metagene.normalize_positions(
    annotated_sites: pl.DataFrame,
    split_strategy: str = "median",
    bin_number: int = 100,
    weight_col_index: list[int] | None = None
) -> tuple[pl.DataFrame, dict, tuple]
```

**Parameters:**
- `annotated_sites`: DataFrame from map_to_transcripts() with transcript coordinates
- `split_strategy`: Strategy for calculating gene region splits ("median" or other)
- `bin_number`: Number of bins for discretizing the normalized positions (default: 100)
- `weight_col_index`: Optional list of column indices for weighting

**Returns:**
- Tuple of (gene_bins, gene_stats, gene_splits):
  - `gene_bins`: DataFrame with normalized position bins and counts
  - `gene_stats`: Dictionary with statistics for each feature type
  - `gene_splits`: Tuple of (5'UTR ratio, CDS ratio, 3'UTR ratio)

**Example:**
```python
gene_bins, gene_stats, gene_splits = metagene.normalize_positions(
    annotated_sites, 
    split_strategy="median", 
    bin_number=100
)
print(f"5'UTR: {gene_splits[0]:.3f}, CDS: {gene_splits[1]:.3f}, 3'UTR: {gene_splits[2]:.3f}")
```

### show_summary_stats()

Display summary statistics for the analysis.

```python
metagene.show_summary_stats(data: PyRanges) -> None
```

**Parameters:**
- `data`: PyRanges with analysis results

---

## Plotting Functions

### plot_profile()

Generate a metagene profile plot.

```python
metagene.plot_profile(
    gene_bins: pl.DataFrame,
    gene_splits: tuple[float, float, float],
    output_path: str,
    figsize: tuple[int, int] = (10, 5)
) -> None
```

**Parameters:**
- `gene_bins`: DataFrame with binned data from normalize_positions()
- `gene_splits`: Tuple of gene region ratios from normalize_positions()
- `output_path`: Path for output image file
- `figsize`: Figure size (width, height) in inches

**Example:**
```python
gene_bins, gene_stats, gene_splits = metagene.normalize_positions(annotated_sites)
metagene.plot_profile(gene_bins, gene_splits, "metagene_plot.png")
```
- `alpha`: Line transparency

### plot_profile()

Create detailed metagene profile with customization options.

```python
metagene.plot_profile(
    data: PyRanges,
    output_path: str,
    **kwargs
) -> None
```

### plot_binned_statistics()

Generate binned statistics plot.

```python
metagene.plot_binned_statistics(
    data: PyRanges,
    output_path: str,
    **kwargs
) -> None
```

---

## Utility Functions

### annotate_with_features()

Annotate sites with overlapping genomic features.

```python
metagene.annotate_with_features(
    sites: PyRanges,
    features: PyRanges
) -> PyRanges
```

### calculate_bin_statistics()

Calculate statistics for binned data.

```python
metagene.calculate_bin_statistics(
    data: PyRanges,
    bins: int = 100
) -> PyRanges
```

---

## Command Line Interface

### Main Command

```bash
metagene [OPTIONS]
```

### Options

| Option | Type | Description |
|--------|------|-------------|
| `-i, --input` | PATH | Input file path |
| `-o, --output` | PATH | Output file path |
| `-r, --reference` | TEXT | Built-in reference (e.g., GRCh38) |
| `-g, --gtf` | PATH | Custom GTF file |
| `-p, --output-figure` | PATH | Output plot file |
| `--region` | CHOICE | Region to analyze (all/5utr/cds/3utr) |
| `--bins` | INTEGER | Number of bins (default: 100) |
| `--with-header` | FLAG | Input file has header |
| `-m, --meta-columns` | TEXT | Column indices for coordinates |
| `--list` | FLAG | List available references |
| `--download` | TEXT | Download reference |

### Examples

```bash
# Basic analysis
metagene -i sites.bed -o results.tsv -r GRCh38

# With custom parameters
metagene -i sites.tsv --with-header -m "1,2,3" -r GRCm39 --bins 200

# List and download references
metagene --list
metagene --download GRCh38
```

---

## Data Structures

### PyRanges Objects

The package uses PyRanges objects to represent genomic intervals and annotations. Key columns include:

- `Chromosome`: Chromosome name
- `Start`: Start position (0-based)
- `End`: End position (exclusive)
- `Strand`: Strand ('+' or '-')

Additional columns may be present depending on the analysis step:

- `transcript_id`: Transcript identifier
- `gene_id`: Gene identifier
- `normalized_position`: Position normalized to gene structure
- `bin`: Bin number for aggregated data

---

## Error Handling

The package provides informative error messages for common issues:

- **File not found**: Clear message with suggested solutions
- **Invalid format**: Description of expected format
- **Missing references**: Automatic download prompts
- **Column specification**: Helpful guidance for meta-columns parameter

---

## Type Hints

All functions include comprehensive type hints for better IDE support and code clarity. Import types as needed:

```python
from typing import List, Optional, Union
import pyranges as pr
```
