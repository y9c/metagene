#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Metagene Analysis Package - Refactored for PyRanges v1

# Core modules with clean API names
from .gtf import prepare_exon_ref, load_gtf
from .io import load_sites, parse_feature_file, load_reference
from .overlap import annotate_with_features, calculate_bin_statistics
from .annotation import map_to_transcripts, normalize_positions, show_summary_stats
from .plotting import plot_profile


# Export main functions
__all__ = [
    # Core analysis functions - new clean API
    "annotate_with_features",
    "calculate_bin_statistics",
    "map_to_transcripts",
    "normalize_positions",
    "show_summary_stats",
    # Data I/O functions - new clean API
    "load_sites",
    "parse_feature_file",
    "load_reference",
    # GTF parsing functions - new clean API
    "prepare_exon_ref",
    "load_gtf",
    # Plotting functions - new clean API
    "plot_profile",
]
