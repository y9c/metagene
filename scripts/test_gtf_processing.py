#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test GTF Processing Script

This script tests the GTF processing functionality with a small example
to ensure everything works before processing large files.
"""

import os
import sys
import tempfile
from pathlib import Path  # noqa: F401  (may be useful later)

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from process_gtf_to_parquet import process_gtf_file


def create_test_gtf():
    """Create a test GTF file for validation."""

    gtf_content = """##gtf-version 2.2
##provider: GENCODE
##format: gtf
##date: 2023-01-01
chr1	GENCODE	gene	1000	5000	.	+	.	gene_id "ENSG00000001"; gene_name "GENE1"; gene_type "protein_coding";
chr1	GENCODE	transcript	1000	5000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1"; transcript_name "GENE1-001"; transcript_type "protein_coding";
chr1	GENCODE	exon	1000	1500	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "1";
chr1	GENCODE	exon	2000	2500	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "2";
chr1	GENCODE	exon	3000	5000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "3";
chr1	GENCODE	start_codon	1200	1202	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001";
chr1	GENCODE	stop_codon	4800	4802	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001";
chr1	GENCODE	gene	10000	15000	.	-	.	gene_id "ENSG00000002"; gene_name "GENE2"; gene_type "protein_coding";
chr1	GENCODE	transcript	10000	15000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; gene_name "GENE2"; transcript_name "GENE2-001"; transcript_type "protein_coding";
chr1	GENCODE	exon	10000	12000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; exon_number "1";
chr1	GENCODE	exon	13000	13500	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; exon_number "2";
chr1	GENCODE	exon	14000	15000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; exon_number "3";
chr1	GENCODE	start_codon	14800	14802	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002";
chr1	GENCODE	stop_codon	10200	10202	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002";
"""

    # Create temporary GTF file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False) as f:
        f.write(gtf_content)
        return f.name


def test_processing():
    """Test the GTF processing pipeline."""

    print("Creating test GTF file...")
    test_gtf = create_test_gtf()

    try:
        print(f"Test GTF file created: {test_gtf}")

        # Process the test file
        output_path = process_gtf_file(test_gtf, "test_output")

        assert os.path.exists(output_path), "Output parquet file should be created"

        # Test loading the file
        import polars as pl

        try:
            df = pl.read_parquet(output_path)
        finally:
            # Clean up output file regardless of read success
            if os.path.exists(output_path):
                os.remove(output_path)
                print(f"\nCleaned up: {output_path}")

        # Basic sanity checks on dataframe
        assert df.height > 0, "Parquet should contain rows"
        assert {"Chromosome", "Start", "End"}.issubset(set(df.columns)), (
            "Expected columns missing"
        )

    finally:
        # Clean up test GTF file
        if os.path.exists(test_gtf):
            os.remove(test_gtf)
            print(f"Cleaned up: {test_gtf}")


if __name__ == "__main__":
    print("=" * 60)
    print("TESTING GTF TO PARQUET CONVERSION")
    print("=" * 60)

    success = test_processing()

    if success:
        print("\n🎉 Test completed successfully!")
        print("The GTF processing pipeline is working correctly.")
        print("You can now process large GTF files with confidence.")
    else:
        print("\n❌ Test failed!")
        print("Please check the error messages above.")
        sys.exit(1)
