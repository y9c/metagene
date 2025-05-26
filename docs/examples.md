---
layout: default
title: Examples
nav_order: 4
---

# Examples
{: .no_toc }

Complete examples and workflows for common use cases.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## RNA Modification Analysis

### m6A Methylation Sites

This example analyzes the distribution of m6A methylation sites across mRNA transcripts.

```python
import metagene

# Load m6A sites from METTL3 knockdown experiment
sites = metagene.load_sites(
    "m6a_sites.bed",
    meta_col_index=[0, 1, 2, 5]  # chr, start, end, strand
)

# Use human reference genome
reference = metagene.load_reference("GRCh38")

# Map sites to transcripts
mapped = metagene.map_to_transcripts(sites, reference)

# Focus on complete transcripts (5'UTR + CDS + 3'UTR)
normalized = metagene.normalize_positions(
    mapped, 
    region="all", 
    bins=300  # Higher resolution
)

# Generate publication-ready plot
metagene.plot_profile(
    normalized,
    "m6a_metagene_profile.png",
    title="m6A Methylation Sites Distribution",
    figsize=(12, 6),
    color="darkred"
)

# Show summary statistics
metagene.show_summary_stats(normalized)
```

### Command Line Version

```bash
# Quick analysis via CLI
metagene -i m6a_sites.bed \
         -o m6a_results.tsv \
         -r GRCh38 \
         -p m6a_plot.png \
         --bins 300 \
         --region all
```

---

## ChIP-seq Peak Analysis

### Transcription Factor Binding

Analyze the distribution of transcription factor binding peaks relative to gene structure.

```python
import metagene

# Load ChIP-seq peaks
peaks = metagene.load_sites(
    "tfbs_peaks.bed",
    meta_col_index=[0, 1, 2, 5]  # BED format
)

# Use mouse reference
reference = metagene.load_reference("GRCm39")

# Map to transcripts
mapped = metagene.map_to_transcripts(peaks, reference)

# Analyze different regions separately
regions = ["5utr", "cds", "3utr"]
colors = ["blue", "green", "orange"]

for region, color in zip(regions, colors):
    normalized = metagene.normalize_positions(
        mapped, 
        region=region, 
        bins=100
    )
    
    metagene.plot_profile(
        normalized,
        f"tfbs_{region}_profile.png",
        title=f"TFBS Distribution - {region.upper()}",
        color=color
    )
```

---

## Multi-Sample Comparison

### Comparing Treatment Conditions

Compare metagene profiles between different experimental conditions.

```python
import metagene
import matplotlib.pyplot as plt

# Load data from different conditions
control_sites = metagene.load_sites("control_sites.tsv", with_header=True)
treated_sites = metagene.load_sites("treated_sites.tsv", with_header=True)

# Use same reference for both
reference = metagene.load_reference("GRCh38")

# Process both datasets
datasets = {
    "Control": control_sites,
    "Treated": treated_sites
}

results = {}
for condition, sites in datasets.items():
    mapped = metagene.map_to_transcripts(sites, reference)
    normalized = metagene.normalize_positions(mapped, region="all")
    results[condition] = normalized

# Create combined plot
fig, ax = plt.subplots(figsize=(12, 6))

colors = ["blue", "red"]
for i, (condition, data) in enumerate(results.items()):
    # Calculate profile (simplified - use actual binning function)
    profile = metagene.calculate_bin_statistics(data)
    ax.plot(profile['bin'], profile['count'], 
            label=condition, color=colors[i], linewidth=2)

ax.set_xlabel("Relative Position")
ax.set_ylabel("Signal Density")
ax.set_title("Metagene Profile Comparison")
ax.legend()
plt.tight_layout()
plt.savefig("comparison_plot.png", dpi=300)
```

---

## Custom GTF Analysis

### Working with Custom Annotations

Use your own GTF file for specialized analysis.

```python
import metagene

# Load sites
sites = metagene.load_sites(
    "custom_sites.bed", 
    meta_col_index=[0, 1, 2, 5]
)

# Load custom GTF (e.g., long non-coding RNAs)
reference = metagene.load_gtf("lncRNA_annotations.gtf")

# Standard analysis workflow
mapped = metagene.map_to_transcripts(sites, reference)
normalized = metagene.normalize_positions(mapped, region="all")

# Generate plot
metagene.plot_profile(
    normalized,
    "lncRNA_analysis.png",
    title="Sites Distribution on lncRNAs"
)
```

---

## High-Resolution Analysis

### Fine-Grained Resolution

For detailed analysis around specific features.

```python
import metagene

# Load high-confidence sites
sites = metagene.load_sites("high_conf_sites.bed")
reference = metagene.load_reference("GRCh38")

# High-resolution analysis
mapped = metagene.map_to_transcripts(sites, reference)

# Focus on CDS with high resolution
normalized = metagene.normalize_positions(
    mapped, 
    region="cds", 
    bins=500  # Very high resolution
)

# Create detailed plot
metagene.plot_profile(
    normalized,
    "high_res_cds.png",
    title="High-Resolution CDS Profile",
    figsize=(15, 8)
)
```

---

## Batch Processing

### Processing Multiple Files

Analyze multiple samples in batch.

```python
import metagene
from pathlib import Path

# Input directory with multiple BED files
input_dir = Path("input_samples/")
output_dir = Path("results/")
output_dir.mkdir(exist_ok=True)

# Load reference once
reference = metagene.load_reference("GRCh38")

# Process each file
for bed_file in input_dir.glob("*.bed"):
    sample_name = bed_file.stem
    
    print(f"Processing {sample_name}...")
    
    # Load and analyze
    sites = metagene.load_sites(str(bed_file))
    mapped = metagene.map_to_transcripts(sites, reference)
    normalized = metagene.normalize_positions(mapped)
    
    # Save results
    output_file = output_dir / f"{sample_name}_results.tsv"
    plot_file = output_dir / f"{sample_name}_plot.png"
    
    # Save normalized data (implementation depends on format)
    # normalized.to_csv(output_file, sep='\t')
    
    metagene.plot_profile(
        normalized,
        str(plot_file),
        title=f"Metagene Profile - {sample_name}"
    )
    
    print(f"Results saved to {output_file}")
```

### Command Line Batch Processing

```bash
#!/bin/bash
# Batch processing script

for file in input_samples/*.bed; do
    basename=$(basename "$file" .bed)
    echo "Processing $basename..."
    
    metagene -i "$file" \
             -o "results/${basename}_results.tsv" \
             -r GRCh38 \
             -p "results/${basename}_plot.png"
done
```

---

## Quality Control

### Data Validation and QC

Perform quality control checks on your analysis.

```python
import metagene

# Load data
sites = metagene.load_sites("sites.bed")
reference = metagene.load_reference("GRCh38")

# Basic statistics
print(f"Total sites: {len(sites)}")
print(f"Chromosomes: {sites.df['Chromosome'].unique()}")

# Map to transcripts
mapped = metagene.map_to_transcripts(sites, reference)

# Check mapping success rate
total_sites = len(sites)
mapped_sites = len(mapped)
mapping_rate = mapped_sites / total_sites * 100

print(f"Mapping rate: {mapping_rate:.1f}%")

# Check strand distribution
strand_counts = mapped.df['Strand'].value_counts()
print(f"Strand distribution: {strand_counts}")

# Proceed with analysis only if mapping rate is acceptable
if mapping_rate > 50:  # Threshold
    normalized = metagene.normalize_positions(mapped)
    metagene.plot_profile(normalized, "qc_passed_plot.png")
else:
    print("Warning: Low mapping rate, check input data quality")
```

---

## Advanced Visualization

### Custom Plot Styling

Create publication-ready plots with custom styling.

```python
import metagene
import matplotlib.pyplot as plt
import seaborn as sns

# Set publication style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Load and process data
sites = metagene.load_sites("sites.bed")
reference = metagene.load_reference("GRCh38")
mapped = metagene.map_to_transcripts(sites, reference)
normalized = metagene.normalize_positions(mapped)

# Create custom plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot with custom styling
# (This would use internal plotting functions)
metagene.plot_profile(
    normalized,
    "publication_ready.png",
    title="Metagene Analysis",
    figsize=(10, 6),
    color="#2E86AB",
    alpha=0.8
)

# Add additional styling
plt.xlabel("Relative Position", fontsize=12, fontweight='bold')
plt.ylabel("Signal Density", fontsize=12, fontweight='bold')
plt.title("RNA Modification Sites Distribution", fontsize=14, fontweight='bold')
plt.tight_layout()
```

---

## Tips and Best Practices

### Memory Optimization

For large datasets, consider these optimization strategies:

```python
import metagene

# For very large files, process in chunks
def process_large_file(filename, chunk_size=10000):
    # This is a conceptual example
    # Actual implementation would depend on file format
    
    reference = metagene.load_reference("GRCh38")
    results = []
    
    # Process file in chunks
    for chunk in read_file_chunks(filename, chunk_size):
        sites = metagene.load_sites(chunk)
        mapped = metagene.map_to_transcripts(sites, reference)
        normalized = metagene.normalize_positions(mapped)
        results.append(normalized)
    
    # Combine results
    combined = combine_results(results)
    return combined
```

### Parameter Optimization

Finding optimal parameters for your analysis:

```python
import metagene

# Test different bin numbers
bin_sizes = [50, 100, 200, 500]
sites = metagene.load_sites("sites.bed")
reference = metagene.load_reference("GRCh38")
mapped = metagene.map_to_transcripts(sites, reference)

for bins in bin_sizes:
    normalized = metagene.normalize_positions(mapped, bins=bins)
    metagene.plot_profile(
        normalized,
        f"test_bins_{bins}.png",
        title=f"Bins: {bins}"
    )
```

These examples demonstrate the flexibility and power of the Metagene package for various genomic analysis scenarios. Adapt them to your specific research needs!
