#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test cases for map_to_local function

Tests the map_to_local functionality that converts global genomic coordinates
to local reference coordinates, using ruranges instead of pyranges1.
"""

import pytest
import polars as pl
import numpy as np
from metagene.map_to_local import map_to_local


class TestMapToLocal:
    """Test suite for map_to_local function"""
    
    def test_basic_mapping_plus_strand(self):
        """Test basic coordinate mapping on plus strand"""
        # Create query intervals
        query = pl.DataFrame({
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [100, 150, 250],
            "End": [110, 160, 260],
            "Strand": ["+", "+", "+"],
        })
        
        # Create reference (single exon)
        reference = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [90],
            "End": [300],
            "Strand": ["+"],
            "transcript_id": ["TX1"],
        })
        
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        
        # Check that coordinates are transformed correctly
        assert result.height == 3
        assert result["Chromosome"].to_list() == ["TX1", "TX1", "TX1"]
        # On plus strand: local_start = global_start - ref_start = 100 - 90 = 10
        assert result["Start"][0] == 10
        assert result["End"][0] == 20
        assert result["Start"][1] == 60
        assert result["End"][1] == 70
        
    def test_basic_mapping_minus_strand(self):
        """Test basic coordinate mapping on minus strand"""
        query = pl.DataFrame({
            "Chromosome": ["chr1", "chr1"],
            "Start": [100, 150],
            "End": [110, 160],
            "Strand": ["-", "-"],
        })
        
        reference = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [90],
            "End": [300],
            "Strand": ["-"],
            "transcript_id": ["TX1"],
        })
        
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        
        # Check that coordinates are transformed correctly
        assert result.height == 2
        assert result["Chromosome"].to_list() == ["TX1", "TX1"]
        # On minus strand: coordinates are flipped
        # local_start = ref_end - query_end = 300 - 110 = 190
        assert result["Start"][0] == 190
        assert result["End"][0] == 200
    
    def test_multi_exon_transcript(self):
        """Test mapping across multiple exons"""
        query = pl.DataFrame({
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [100, 250, 450],
            "End": [110, 260, 460],
            "Strand": ["+", "+", "+"],
            "query_id": ["Q1", "Q2", "Q3"],  # Add ID to track queries
        })
        
        # Two exons in the transcript
        reference = pl.DataFrame({
            "Chromosome": ["chr1", "chr1"],
            "Start": [90, 240],
            "End": [200, 500],
            "Strand": ["+", "+"],
            "transcript_id": ["TX1", "TX1"],
        })
        
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        
        # All queries should map to TX1
        assert result.height == 3
        assert all(result["Chromosome"] == "TX1")
        
        # Sort by query_id to get consistent order
        result = result.sort("query_id")
        
        # First query (Q1) in first exon: 100 - 90 = 10
        q1 = result.filter(pl.col("query_id") == "Q1")
        assert q1["Start"][0] == 10
        assert q1["End"][0] == 20
        
        # Second query (Q2) in second exon: cumsum_start = 110, local = 110 + (250 - 240) = 120
        q2 = result.filter(pl.col("query_id") == "Q2")
        assert q2["Start"][0] == 120
        assert q2["End"][0] == 130
        
        # Third query (Q3) in second exon: cumsum_start = 110, local = 110 + (450 - 240) = 320
        q3 = result.filter(pl.col("query_id") == "Q3")
        assert q3["Start"][0] == 320
        assert q3["End"][0] == 330
    
    def test_strand_transformation(self):
        """Test that strand is transformed correctly"""
        # Query on plus strand, reference on minus strand
        query = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [110],
            "Strand": ["+"],
        })
        
        reference = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [90],
            "End": [300],
            "Strand": ["-"],
            "transcript_id": ["TX1"],
        })
        
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        
        # When query and reference have different strands, result should be "-"
        assert result["Strand"][0] == "-"
        
    def test_no_overlap(self):
        """Test behavior when there's no overlap"""
        query = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [1000],
            "End": [1100],
            "Strand": ["+"],
        })
        
        reference = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [90],
            "End": [300],
            "Strand": ["+"],
            "transcript_id": ["TX1"],
        })
        
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        
        # Should return empty dataframe
        assert result.height == 0
    
    def test_keep_global_coordinates(self):
        """Test keeping global coordinates in result"""
        query = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [110],
            "Strand": ["+"],
        })
        
        reference = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [90],
            "End": [300],
            "Strand": ["+"],
            "transcript_id": ["TX1"],
        })
        
        result = map_to_local(
            query, 
            reference, 
            ref_id_col="transcript_id",
            keep_global_chrom=True,
            keep_global_loc=True
        )
        
        # Check that global coordinates are preserved
        assert "Chromosome_global" in result.columns
        assert "Start_global" in result.columns
        assert "End_global" in result.columns
        assert result["Chromosome_global"][0] == "chr1"
        assert result["Start_global"][0] == 100
        assert result["End_global"][0] == 110
    
    def test_match_by_filter(self):
        """Test filtering by match_by columns"""
        query = pl.DataFrame({
            "Chromosome": ["chr1", "chr1"],
            "Start": [100, 150],
            "End": [110, 160],
            "Strand": ["+", "+"],
            "gene_id": ["G1", "G2"],
        })
        
        reference = pl.DataFrame({
            "Chromosome": ["chr1", "chr1"],
            "Start": [90, 140],
            "End": [200, 250],
            "Strand": ["+", "+"],
            "transcript_id": ["TX1", "TX2"],
            "gene_id": ["G1", "G2"],
        })
        
        result = map_to_local(
            query, 
            reference, 
            ref_id_col="transcript_id",
            match_by="gene_id"
        )
        
        # Should only match queries with their corresponding gene_id
        assert result.height == 2
        # First query (G1) should map to TX1
        # Second query (G2) should map to TX2
    
    def test_multiple_references_same_chromosome(self):
        """Test when multiple references exist on same chromosome"""
        query = pl.DataFrame({
            "Chromosome": ["chr1", "chr1"],
            "Start": [100, 300],
            "End": [110, 310],
            "Strand": ["+", "+"],
        })
        
        # Two separate transcripts on chr1
        reference = pl.DataFrame({
            "Chromosome": ["chr1", "chr1"],
            "Start": [90, 290],
            "End": [200, 400],
            "Strand": ["+", "+"],
            "transcript_id": ["TX1", "TX2"],
        })
        
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        
        # First query should map to TX1, second to TX2
        assert result.height == 2
        assert result["Chromosome"][0] == "TX1"
        assert result["Chromosome"][1] == "TX2"
    
    def test_additional_query_columns_preserved(self):
        """Test that additional columns in query are preserved"""
        query = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [110],
            "Strand": ["+"],
            "score": [0.95],
            "name": ["SNP1"],
        })
        
        reference = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [90],
            "End": [300],
            "Strand": ["+"],
            "transcript_id": ["TX1"],
        })
        
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        
        # Additional columns should be preserved
        assert "score" in result.columns
        assert "name" in result.columns
        assert result["score"][0] == 0.95
        assert result["name"][0] == "SNP1"
    
    def test_missing_required_columns_error(self):
        """Test that missing required columns raise appropriate errors"""
        # Query missing Strand column
        query = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [110],
        })
        
        reference = pl.DataFrame({
            "Chromosome": ["chr1"],
            "Start": [90],
            "End": [300],
            "Strand": ["+"],
            "transcript_id": ["TX1"],
        })
        
        with pytest.raises(ValueError, match="Query DataFrame must have 'Strand' column"):
            map_to_local(query, reference, ref_id_col="transcript_id")
    
    def test_real_world_example_snp_mapping(self):
        """Test with a realistic SNP mapping example"""
        # SNPs to map
        snps = pl.DataFrame({
            "Chromosome": ["chr1", "chr1", "chr1", "chr2"],
            "Start": [1000, 1500, 2000, 5000],
            "End": [1001, 1501, 2001, 5001],
            "Strand": ["+", "+", "+", "+"],
            "rsid": ["rs1", "rs2", "rs3", "rs4"],
        })
        
        # Transcript annotations with exons
        transcripts = pl.DataFrame({
            "Chromosome": ["chr1", "chr1", "chr2"],
            "Start": [900, 1900, 4900],
            "End": [1600, 2200, 5200],
            "Strand": ["+", "+", "+"],
            "transcript_id": ["ENST001", "ENST001", "ENST002"],
        })
        
        result = map_to_local(snps, transcripts, ref_id_col="transcript_id")
        
        # Should map SNPs that overlap transcripts
        assert result.height >= 3  # At least rs1, rs2, rs3 should map
        assert all(result["Chromosome"].is_in(["ENST001", "ENST002"]))
        
        # Check that rsid is preserved
        assert "rsid" in result.columns


def test_compare_with_demo_data():
    """
    Test map_to_local with actual demo data from the repository
    to ensure it produces sensible results
    """
    from pathlib import Path
    from metagene import load_gtf, load_sites
    
    # Load demo data
    test_dir = Path(__file__).parent
    gtf_file = test_dir / "example.gtf.gz"
    sites_file = test_dir / "sites.tsv.gz"
    
    if not gtf_file.exists() or not sites_file.exists():
        pytest.skip("Demo data files not found")
    
    # Load exon reference
    exon_ref = load_gtf(str(gtf_file))
    
    # Load input sites
    input_sites = load_sites(
        str(sites_file), 
        with_header=True, 
        meta_col_index=[0, 1, 2]
    )
    
    # Map to local coordinates using map_to_local
    result = map_to_local(
        input_sites, 
        exon_ref, 
        ref_id_col="transcript_id",
        keep_global_loc=True
    )
    
    # Basic sanity checks
    assert result.height > 0, "Should have some mapped results"
    assert "Chromosome" in result.columns
    assert "Start" in result.columns
    assert "End" in result.columns
    
    # All Start should be >= 0 (local coordinates)
    assert all(result["Start"] >= 0), "Local start coordinates should be non-negative"
    
    # End should be > Start
    assert all(result["End"] > result["Start"]), "End should be greater than Start"
    
    # If we kept global locations, they should be preserved
    if "Start_global" in result.columns:
        assert all(result["Start_global"] >= 0)
        assert all(result["End_global"] > result["Start_global"])
    
    print(f"\nMapped {result.height} sites to local transcript coordinates")
    print(f"Unique transcripts: {result['Chromosome'].n_unique()}")
    print(f"Local coordinate range: {result['Start'].min()} to {result['End'].max()}")


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "-s"])
