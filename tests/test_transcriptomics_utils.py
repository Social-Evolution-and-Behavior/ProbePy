"""
Tests for probepy.transcriptomics.utils module.

This module tests utility functions for transcriptome processing including
sequence extraction and GTF file parsing.
"""

import pytest
import pandas as pd
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, MagicMock, patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from probepy.transcriptomics.utils import get_sequence, gtf_to_dataframe


class TestGetSequence:
    """Test get_sequence function for extracting genomic sequences."""
    
    @pytest.fixture
    def mock_genome_dict(self):
        """Create a mock genome sequence dictionary."""
        # Create mock chromosome sequences
        chr1_seq = Seq("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
        chr2_seq = Seq("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")
        
        seq_dict = {
            "chr1": SeqRecord(chr1_seq, id="chr1"),
            "chr2": SeqRecord(chr2_seq, id="chr2")
        }
        
        return seq_dict
    
    def test_basic_forward_strand(self, mock_genome_dict):
        """Test basic sequence extraction on forward strand."""
        # Extract first 10 bases from chr1
        result = get_sequence(mock_genome_dict, "chr1", 0, 10, "+")
        
        assert result == "ATCGATCGAT"
        assert isinstance(result, str)
        assert len(result) == 10
    
    def test_basic_reverse_strand(self, mock_genome_dict):
        """Test sequence extraction on reverse strand."""
        # Extract and reverse complement
        result = get_sequence(mock_genome_dict, "chr1", 0, 10, "-")
        
        # Should be reverse complement of ATCGATCGAT = ATCGATCGAT
        assert len(result) == 10
        assert isinstance(result, str)
        # Verify it's different from forward (except for palindromes)
        forward = str(mock_genome_dict["chr1"].seq[0:10])
        assert result == str(Seq(forward).reverse_complement())
    
    def test_multiple_chromosomes(self, mock_genome_dict):
        """Test extraction from different chromosomes."""
        seq1 = get_sequence(mock_genome_dict, "chr1", 0, 5, "+")
        seq2 = get_sequence(mock_genome_dict, "chr2", 0, 5, "+")
        
        assert seq1 == "ATCGA"
        assert seq2 == "GCGCG"
        assert seq1 != seq2
    
    def test_different_positions(self, mock_genome_dict):
        """Test extraction from different positions."""
        seq1 = get_sequence(mock_genome_dict, "chr1", 0, 5, "+")
        seq2 = get_sequence(mock_genome_dict, "chr1", 5, 10, "+")
        seq3 = get_sequence(mock_genome_dict, "chr1", 10, 15, "+")
        
        assert seq1 == "ATCGA"
        assert seq2 == "TCGAT"
        assert seq3 == "CGATC"
        # Verify they're consecutive
        assert seq1 + seq2 + seq3 == "ATCGATCGATCGATC"
    
    def test_uppercase_output(self, mock_genome_dict):
        """Test that output is uppercase."""
        result = get_sequence(mock_genome_dict, "chr1", 0, 10, "+")
        
        assert result.isupper()
        assert result == result.upper()
    
    def test_invalid_chromosome_raises_error(self, mock_genome_dict):
        """Test that invalid chromosome raises ValueError."""
        with pytest.raises(ValueError, match="Chromosome .* not found"):
            get_sequence(mock_genome_dict, "chr99", 0, 10, "+")
    
    def test_start_greater_than_end_raises_error(self, mock_genome_dict):
        """Test that start >= end raises ValueError."""
        with pytest.raises(ValueError, match="Start position .* is greater than or equal to end"):
            get_sequence(mock_genome_dict, "chr1", 10, 10, "+")
        
        with pytest.raises(ValueError, match="Start position .* is greater than or equal to end"):
            get_sequence(mock_genome_dict, "chr1", 15, 10, "+")
    
    def test_edge_cases_sequence_boundaries(self, mock_genome_dict):
        """Test extraction at sequence boundaries."""
        chr_length = len(mock_genome_dict["chr1"].seq)
        
        # Extract from end of chromosome
        result = get_sequence(mock_genome_dict, "chr1", chr_length - 5, chr_length, "+")
        assert len(result) == 5
        
        # Extract entire chromosome
        result = get_sequence(mock_genome_dict, "chr1", 0, chr_length, "+")
        assert len(result) == chr_length
    
    def test_reverse_complement_correctness(self, mock_genome_dict):
        """Test that reverse complement is calculated correctly."""
        # Extract known sequence
        forward = get_sequence(mock_genome_dict, "chr1", 0, 4, "+")
        reverse = get_sequence(mock_genome_dict, "chr1", 0, 4, "-")
        
        # Manual reverse complement of ATCG should be CGAT
        assert forward == "ATCG"
        assert reverse == "CGAT"
    
    def test_single_base_extraction(self, mock_genome_dict):
        """Test extraction of single base."""
        result = get_sequence(mock_genome_dict, "chr1", 0, 1, "+")
        
        assert len(result) == 1
        assert result == "A"
    
    def test_long_sequence_extraction(self, mock_genome_dict):
        """Test extraction of longer sequences."""
        result = get_sequence(mock_genome_dict, "chr1", 0, 48, "+")
        
        assert len(result) == 48
        assert isinstance(result, str)


class TestGtfToDataframe:
    """Test gtf_to_dataframe function for parsing GTF files."""
    
    @pytest.fixture
    def simple_gtf_file(self):
        """Create a simple GTF file for testing."""
        gtf_content = """# Sample GTF file
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene1";
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene1"; transcript_id "TRANS001";
chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene1"; transcript_id "TRANS001";
chr1\tHAVANA\tCDS\t1020\t1080\t.\t+\t0\tgene_id "GENE001"; gene_name "TestGene1"; transcript_id "TRANS001";
chr2\tENSEMBL\tgene\t3000\t4000\t.\t-\t.\tgene_id "GENE002"; gene_name "TestGene2";
chr2\tENSEMBL\ttranscript\t3000\t4000\t.\t-\t.\tgene_id "GENE002"; gene_name "TestGene2"; transcript_id "TRANS002";
chr2\tENSEMBL\texon\t3000\t3500\t.\t-\t.\tgene_id "GENE002"; gene_name "TestGene2"; transcript_id "TRANS002";
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            temp_path = f.name
        
        yield temp_path
        
        # Cleanup
        try:
            os.unlink(temp_path)
        except (OSError, FileNotFoundError):
            pass
    
    @pytest.fixture
    def gtf_with_comments(self):
        """Create GTF file with comments and empty lines."""
        gtf_content = """##gff-version 3
# This is a comment
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene1";

# Another comment
chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            temp_path = f.name
        
        yield temp_path
        
        try:
            os.unlink(temp_path)
        except (OSError, FileNotFoundError):
            pass
    
    def test_basic_parsing(self, simple_gtf_file):
        """Test basic GTF file parsing."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # Check DataFrame structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        
        # Check expected columns
        expected_columns = ['seqname', 'source', 'feature', 'start', 'end', 
                          'score', 'strand', 'frame', 'gene_id', 'gene_name', 
                          'transcript_id', 'attributes']
        for col in expected_columns:
            assert col in df.columns
    
    def test_data_types(self, simple_gtf_file):
        """Test that columns have correct data types."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # start and end should be integers
        assert df['start'].dtype == int
        assert df['end'].dtype == int
        
        # Other columns should be objects/strings
        assert df['seqname'].dtype == object
        assert df['feature'].dtype == object
        assert df['strand'].dtype == object
    
    def test_gene_id_extraction(self, simple_gtf_file):
        """Test gene_id extraction from attributes."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # Check that gene_id is extracted correctly
        assert 'gene_id' in df.columns
        assert 'GENE001' in df['gene_id'].values
        assert 'GENE002' in df['gene_id'].values
    
    def test_transcript_id_extraction(self, simple_gtf_file):
        """Test transcript_id extraction from attributes."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # Check that transcript_id is extracted
        assert 'transcript_id' in df.columns
        
        # Transcript entries should have transcript_id
        transcripts = df[df['feature'] == 'transcript']
        assert len(transcripts) > 0
        assert transcripts['transcript_id'].notna().all()
    
    def test_gene_name_extraction(self, simple_gtf_file):
        """Test gene_name extraction from attributes."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        assert 'gene_name' in df.columns
        assert 'TestGene1' in df['gene_name'].values
        assert 'TestGene2' in df['gene_name'].values
    
    def test_coordinates_parsing(self, simple_gtf_file):
        """Test that coordinates are parsed correctly."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # Check gene coordinates
        gene1 = df[(df['feature'] == 'gene') & (df['gene_id'] == 'GENE001')]
        assert len(gene1) > 0
        assert gene1.iloc[0]['start'] == 1000
        assert gene1.iloc[0]['end'] == 2000
        assert gene1.iloc[0]['strand'] == '+'
        
        gene2 = df[(df['feature'] == 'gene') & (df['gene_id'] == 'GENE002')]
        assert len(gene2) > 0
        assert gene2.iloc[0]['start'] == 3000
        assert gene2.iloc[0]['end'] == 4000
        assert gene2.iloc[0]['strand'] == '-'
    
    def test_feature_types(self, simple_gtf_file):
        """Test that different feature types are captured."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        features = df['feature'].unique()
        assert 'gene' in features
        assert 'transcript' in features
        assert 'exon' in features
        assert 'CDS' in features
    
    def test_comment_filtering(self, gtf_with_comments):
        """Test that comments are filtered out."""
        df = gtf_to_dataframe(gtf_with_comments)
        
        # Should only have data rows, not comments
        assert len(df) == 2  # Only 2 non-comment lines
        assert df['feature'].isin(['gene', 'exon']).all()
    
    def test_attributes_dictionary(self, simple_gtf_file):
        """Test that attributes are stored as dictionary."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # Check that attributes column contains dictionaries
        assert 'attributes' in df.columns
        
        first_row_attrs = df.iloc[0]['attributes']
        assert isinstance(first_row_attrs, dict)
        assert 'gene_id' in first_row_attrs
        assert 'gene_name' in first_row_attrs
    
    def test_missing_gene_id_raises_error(self):
        """Test that missing gene_id raises ValueError."""
        # Create GTF file without gene_id
        gtf_content = """chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_name "TestGene";
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            temp_path = f.name
        
        try:
            with pytest.raises(ValueError, match="Some entries are missing gene_id"):
                gtf_to_dataframe(temp_path)
        finally:
            try:
                os.unlink(temp_path)
            except:
                pass
    
    def test_multiple_genes_and_transcripts(self, simple_gtf_file):
        """Test parsing file with multiple genes and transcripts."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # Should have 2 unique genes
        unique_genes = df['gene_id'].unique()
        assert len(unique_genes) == 2
        
        # Should have 2 transcripts
        transcripts = df[df['feature'] == 'transcript']
        assert len(transcripts) == 2
    
    def test_chromosome_names(self, simple_gtf_file):
        """Test that chromosome names are parsed correctly."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        chromosomes = df['seqname'].unique()
        assert 'chr1' in chromosomes
        assert 'chr2' in chromosomes
    
    def test_strand_information(self, simple_gtf_file):
        """Test that strand information is captured."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        # Should have both + and - strands
        strands = df['strand'].unique()
        assert '+' in strands
        assert '-' in strands
    
    def test_source_field(self, simple_gtf_file):
        """Test that source field is parsed."""
        df = gtf_to_dataframe(simple_gtf_file)
        
        sources = df['source'].unique()
        assert 'HAVANA' in sources
        assert 'ENSEMBL' in sources
    
    def test_empty_gtf_file(self):
        """Test behavior with empty GTF file."""
        gtf_content = """# Only comments
# No data
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            temp_path = f.name
        
        try:
            df = gtf_to_dataframe(temp_path)
            # Should return empty DataFrame with correct columns
            assert len(df) == 0
            assert 'gene_id' in df.columns
        finally:
            try:
                os.unlink(temp_path)
            except:
                pass
    
    def test_malformed_attributes(self):
        """Test handling of malformed attribute strings."""
        gtf_content = """chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "GENE001" gene_name "TestGene";
chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tgene_id "GENE001"; transcript_id;
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            temp_path = f.name
        
        try:
            # Should handle gracefully or raise appropriate error
            df = gtf_to_dataframe(temp_path)
            # At minimum should parse the valid entries
            assert 'GENE001' in df['gene_id'].values
        finally:
            try:
                os.unlink(temp_path)
            except:
                pass
    
    def test_utr_features(self):
        """Test parsing of UTR features."""
        gtf_content = """chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";
chr1\tHAVANA\t5UTR\t1000\t1050\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";
chr1\tHAVANA\t3UTR\t1950\t2000\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            temp_path = f.name
        
        try:
            df = gtf_to_dataframe(temp_path)
            features = df['feature'].unique()
            assert '5UTR' in features or 'five_prime_utr' in features
            assert '3UTR' in features or 'three_prime_utr' in features
        finally:
            try:
                os.unlink(temp_path)
            except:
                pass


class TestIntegration:
    """Integration tests for transcriptomics utils."""
    
    def test_gtf_and_sequence_extraction_workflow(self):
        """Test workflow of parsing GTF and extracting sequences."""
        # Create test data
        gtf_content = """chr1\tHAVANA\tgene\t10\t30\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";
chr1\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            gtf_path = f.name
        
        try:
            # Parse GTF
            df = gtf_to_dataframe(gtf_path)
            
            # Create mock genome
            genome_seq = Seq("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
            seq_dict = {"chr1": SeqRecord(genome_seq, id="chr1")}
            
            # Extract sequence for exon
            exon = df[df['feature'] == 'exon'].iloc[0]
            sequence = get_sequence(
                seq_dict,
                exon['seqname'],
                exon['start'] - 1,  # Convert to 0-based
                exon['end'],
                exon['strand']
            )
            
            # Verify
            assert len(sequence) == 10
            assert isinstance(sequence, str)
            
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass
    
    def test_multi_exon_gene_extraction(self):
        """Test extracting sequences for multi-exon gene."""
        gtf_content = """chr1\tHAVANA\tgene\t10\t40\t.\t+\t.\tgene_id "GENE001"; gene_name "TestGene";
chr1\tHAVANA\texon\t10\t15\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";
chr1\tHAVANA\texon\t20\t25\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";
chr1\tHAVANA\texon\t30\t40\t.\t+\t.\tgene_id "GENE001"; transcript_id "TRANS001";
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write(gtf_content)
            gtf_path = f.name
        
        try:
            # Parse GTF
            df = gtf_to_dataframe(gtf_path)
            
            # Create mock genome
            genome_seq = Seq("A" * 50)
            seq_dict = {"chr1": SeqRecord(genome_seq, id="chr1")}
            
            # Extract all exon sequences
            exons = df[df['feature'] == 'exon']
            sequences = []
            for _, exon in exons.iterrows():
                seq = get_sequence(
                    seq_dict,
                    exon['seqname'],
                    exon['start'] - 1,
                    exon['end'],
                    exon['strand']
                )
                sequences.append(seq)
            
            # Verify we got 3 exons
            assert len(sequences) == 3
            assert all(isinstance(s, str) for s in sequences)
            
        finally:
            try:
                os.unlink(gtf_path)
            except:
                pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
