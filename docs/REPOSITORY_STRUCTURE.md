# Metagene Repository - Clean Production Structure

## 📁 Repository Overview

The metagene repository is now production-ready with version 0.0.1, featuring a clean structure, automatic reference downloading, and comprehensive CLI interface.

## 🏗️ Directory Structure

```
metagene/
├── README.md                    # Comprehensive project documentation
├── pyproject.toml              # Project configuration (v0.0.1)
├── Makefile                    # Build and development commands
├── uv.lock                     # Dependency lock file (gitignored)
│
├── .github/                    # GitHub Actions
│   └── workflows/
│       └── publish.yml         # PyPI publishing workflow
│
├── metagene/                   # Main package
│   ├── __init__.py            # Clean API exports
│   ├── gtf.py                 # GTF file processing (PyRanges v1)
│   ├── io.py                  # I/O operations with auto-download
│   ├── overlap.py             # Genomic overlap analysis
│   ├── annotation.py          # Transcript annotation utilities
│   ├── plotting.py            # Visualization functions
│   ├── cli.py                 # Rich CLI with progress bars
│   ├── download.py            # Reference download system
│   ├── config.py              # Configuration and built-in references
│   └── utils.py               # Utility functions
│
├── test/                      # Test files and example data
│   ├── __init__.py
│   ├── test_demo.py           # Complete workflow demo
│   ├── test_builtin.py        # Built-in reference demo
│   ├── example.gtf.gz         # Test GTF file
│   ├── sites.tsv.gz          # Test sites data
│   ├── sites.bed             # Test BED data
│   └── output_*.{tsv,png}    # Test outputs (gitignored)
│
├── scripts/                   # Development scripts
│   ├── README.md             # Script documentation
│   ├── process_gtf_to_parquet.py # GTF processing utility
│   ├── species_mapping_template.py # Species mapping
│   └── test_gtf_processing.py # GTF processing tests
│
└── docs/                     # Documentation
    ├── REPOSITORY_STRUCTURE.md # This file
    ├── metagene_demo_*.png    # Example plots
    └── fig_metagene.svg      # Vector example
```
│
└── legacy/                    # Deprecated code (for reference)
    ├── README.md              # Legacy code documentation
    ├── REFACTORING_SUMMARY.md # Refactoring history
    ├── *.py                   # Old module files
    └── data/                  # Legacy test data
```

## ✅ Clean Repository Features

### 1. **Organized Structure**
- Clear separation between current code (`metagene/`) and legacy code (`legacy/`)
- Documentation centralized in `docs/` folder
- Test files and example data in `test/` folder

### 2. **PyRanges v1 Ready**
- All modules compatible with PyRanges v1.x
- No deprecated `nb_cpu` parameters
- Updated column naming conventions (`_b` suffix)

### 3. **Single-Word Module Names**
- Clean, intuitive module structure
- Easy to import and understand
- Professional naming conventions

### 4. **Backward Compatibility**
- Legacy function names still work via aliases
- Existing code continues to function
- Smooth migration path for users

### 5. **Clean Codebase**
- No Python cache files (`__pycache__/`)
- No temporary files or build artifacts
- No unnecessary documentation in root folder

## 🚀 Usage

### Modern Import Style
```python
from metagene import gtf, io, overlap, analysis, plotting
from metagene.gtf import prepare_exon_ref, read_gtf_with_caching
from metagene.analysis import run_metagene_analysis
```

### Legacy Compatible
```python
from metagene import annotate_with_feature, parse_input, parse_features
from metagene import gtf_to_bed_with_caching  # Legacy function name
```

## 🎯 Repository Status

- ✅ **Clean and organized**
- ✅ **PyRanges v1 compatible**
- ✅ **Fully tested and working**
- ✅ **Documentation complete**
- ✅ **Ready for production use**

This repository structure provides a solid foundation for continued development and maintenance of the metagene package.
