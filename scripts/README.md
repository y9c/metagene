# GTF Processing Scripts for Metagene Built-in References

This directory contains scripts to process large GTF files and convert them to highly compressed Parquet files for built-in species support in the metagene package.

## Overview

The metagene package supports built-in reference annotations for common species. Due to the large size of GTF files (often hundreds of MB to GB), we convert them to highly compressed Parquet files that:

- Reduce file size by 10-50x using ZSTD compression (level 22)
- Enable fast loading with optimized data types
- Maintain all necessary genomic annotation information

## Scripts

### 1. `process_gtf_to_parquet.py`

Main script to convert GTF files to compressed Parquet format. This script now includes all predefined species configurations.

**Usage:**
```bash
# Process all predefined GTF files
python process_gtf_to_parquet.py --all

# Process a single GTF file
python process_gtf_to_parquet.py --single input.gtf output_name

# List all files to be processed
python process_gtf_to_parquet.py --list

# Create species mapping template
python process_gtf_to_parquet.py --template
```

**Examples:**
```bash
# Process all 23 predefined GTF files
python process_gtf_to_parquet.py --all

# Process human GRCh38 GTF
python process_gtf_to_parquet.py --single Homo_sapiens/raw/GRCh38.release110.gtf.gz GRCh38

# List all files that will be processed
python process_gtf_to_parquet.py --list
```

**Features:**
- Includes 23 predefined species and genome builds
- Automatic data type optimization for maximum compression
- Progress reporting and timing information
- Memory usage optimization
- Error handling and validation
- Generates species mapping code for easy integration

### 2. `batch_process_gtf.py` (DEPRECATED)

This script is now deprecated. All functionality has been moved to `process_gtf_to_parquet.py --all`.

### 3. `test_gtf_processing.py`

Test script to validate the processing pipeline.

**Usage:**
```bash
python test_gtf_processing.py
```

**Features:**
- Creates test GTF file
- Validates processing pipeline
- Tests compression and loading
- Cleanup and error reporting

## Common Species Mapping

The following species are commonly supported:

| Species | Output Name | Common Aliases | Description |
|---------|-------------|----------------|-------------|
| Human | `GRCh38` | `hg38` | Human genome GRCh38 |
| Human | `GRCh37` | `hg19` | Human genome GRCh37 |
| Human | `GRCh36` | `hg18` | Human genome GRCh36 |
| Mouse | `GRCm39` | `mm39` | Mouse genome GRCm39 |
| Mouse | `GRCm38` | `mm10` | Mouse genome GRCm38 |
| Mouse | `GRCm37` | `mm9` | Mouse genome GRCm37 |
| Arabidopsis | `TAIR10` | - | Arabidopsis thaliana TAIR10 |
| Rice | `IRGSP-1.0` | - | Rice IRGSP-1.0 |

## Workflow

### On Server with Large GTF Files

1. **Setup environment:**
   ```bash
   # Install required packages
   pip install polars pandas pyranges
   
   # Copy scripts to server
   scp scripts/*.py server:/path/to/gtf/files/
   ```

2. **Test the pipeline:**
   ```bash
   python test_gtf_processing.py
   ```

3. **Process individual files:**
   ```bash
   python process_gtf_to_parquet.py Homo_sapiens.GRCh38.109.gtf GRCh38
   python process_gtf_to_parquet.py Mus_musculus.GRCm39.109.gtf GRCm39
   ```

4. **Or use batch processing:**
   ```bash
   # Create config file
   python batch_process_gtf.py --create-config
   
   # Edit gtf_batch_config.yaml with your file paths
   # Then process all files
   python batch_process_gtf.py
   ```

### Back to Development Machine

1. **Copy Parquet files:**
   ```bash
   scp server:/path/to/gtf/files/*.bed.parquet metagene/data/
   ```

2. **Update metagene package:**
   
   Add to `metagene/io.py`:
   ```python
   BUILTIN_SPECIES = {
       'GRCh38': 'GRCh38.bed.parquet',
       'hg38': 'GRCh38.bed.parquet',
       'GRCm39': 'GRCm39.bed.parquet', 
       'mm39': 'GRCm39.bed.parquet',
       # ... add other species
   }
   ```

3. **Test built-in features:**
   ```python
   from metagene.io import load_builtin_features
   
   # Load human GRCh38 features
   features = load_builtin_features('GRCh38')
   print(f"Loaded {len(features)} features")
   ```

## File Size Examples

Typical compression results:

| Original GTF | Compressed Parquet | Compression Ratio |
|--------------|-------------------|-------------------|
| 1.5 GB | 85 MB | 18x |
| 800 MB | 45 MB | 18x |
| 500 MB | 25 MB | 20x |

## Configuration File Template

```yaml
gtf_files:
  - name: GRCh38
    file: Homo_sapiens.GRCh38.109.gtf
    description: Human genome GRCh38/hg38
  - name: GRCh37
    file: Homo_sapiens.GRCh37.87.gtf
    description: Human genome GRCh37/hg19
  - name: GRCm39
    file: Mus_musculus.GRCm39.109.gtf
    description: Mouse genome GRCm39/mm39

output_directory: ./parquet_files
compression_level: 22
```

## Dependencies

Required Python packages:
- `polars` (for fast data processing and compression)
- `pandas` (for compatibility)
- `pyranges` (for genomic data handling)
- `pyyaml` (for configuration files)

## Troubleshooting

### Memory Issues
- Large GTF files may require substantial RAM
- Consider processing on a high-memory server
- Monitor memory usage during processing

### File Format Issues
- Ensure GTF files follow standard format
- Check for proper gene_id and transcript_id attributes
- Validate exon/start_codon/stop_codon features

### Compression Issues
- ZSTD level 22 provides maximum compression but is slow
- Reduce compression level if speed is more important
- Consider using level 15-18 for good compression/speed balance

## Support

For issues or questions:
1. Run the test script first: `python test_gtf_processing.py`
2. Check error messages for file format issues
3. Verify GTF file integrity and format
4. Ensure sufficient disk space and memory

## Integration with Metagene Package

After generating Parquet files, integrate them into the metagene package:

1. **Copy files to data directory:**
   ```bash
   cp *.bed.parquet /path/to/metagene/data/
   ```

2. **Update io.py with species mapping:**
   ```python
   def load_builtin_features(species: str) -> pd.DataFrame:
       BUILTIN_SPECIES = {
           'GRCh38': 'GRCh38.bed.parquet',
           'hg38': 'GRCh38.bed.parquet',
           # ... other mappings
       }
       
       if species not in BUILTIN_SPECIES:
           available = list(BUILTIN_SPECIES.keys())
           raise ValueError(f"Species '{species}' not available. Available: {available}")
       
       data_dir = os.path.join(os.path.dirname(__file__), "data")
       feature_file = os.path.join(data_dir, BUILTIN_SPECIES[species])
       
       return parse_feature_file(feature_file)
   ```

3. **Test integration:**
   ```python
   import metagene
   features = metagene.load_builtin_features('GRCh38')
   ```
