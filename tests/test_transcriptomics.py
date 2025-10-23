"""
Tests for hcrfish.transcriptomics module.

This module tests the transcriptomics data structures and functions
including Transcriptome, Gene, Transcript, and genomic feature classes.
"""

import pytest
import numpy as np
from typing import Dict, List

from hcrfish.transcriptomics.classes import (
    Transcriptome, Gene, Transcript, Exon, CDS, UTR, Intron
)


class TestTranscriptome:
    """Test Transcriptome class functionality."""
    
    def test_initialization(self):
        """Test transcriptome initialization."""
        transcriptome = Transcriptome()
        
        assert isinstance(transcriptome.genes, dict)
        assert isinstance(transcriptome.transcripts_to_genes, dict)
        assert len(transcriptome.genes) == 0
        assert len(transcriptome.transcripts_to_genes) == 0
    
    def test_add_gene(self, sample_gene):
        """Test adding genes to transcriptome."""
        transcriptome = Transcriptome()
        
        # Add gene
        transcriptome.add_gene(sample_gene)
        
        assert len(transcriptome.genes) == 1
        assert sample_gene.name in transcriptome.genes
        assert transcriptome.genes[sample_gene.name] == sample_gene
    
    def test_get_gene(self, sample_transcriptome):
        """Test gene retrieval."""
        # Get existing gene
        gene = sample_transcriptome.get_gene("Or9a")
        assert gene is not None
        assert gene.name == "Or9a"
        
        # Get non-existing gene
        missing_gene = sample_transcriptome.get_gene("NonExistent")
        assert missing_gene is None
    
    def test_remove_gene(self, sample_transcriptome):
        """Test gene removal."""
        # Remove existing gene
        success = sample_transcriptome.remove_gene("Or9a")
        assert success is True
        assert sample_transcriptome.get_gene("Or9a") is None
        
        # Try to remove non-existing gene
        success = sample_transcriptome.remove_gene("NonExistent")
        assert success is False
    
    def test_search_genes(self, sample_transcriptome):
        """Test gene search functionality."""
        # Search for partial match
        results = sample_transcriptome.search_genes("Or")
        assert len(results) == 1
        assert results[0].name == "Or9a"
        
        # Search case insensitive
        results = sample_transcriptome.search_genes("or")
        assert len(results) == 1
        assert results[0].name == "Or9a"
        
        # Search for multiple matches
        results = sample_transcriptome.search_genes("s")  # Should match "dsx"
        assert len(results) == 1
        assert results[0].name == "dsx"
        
        # Search with no matches
        results = sample_transcriptome.search_genes("xyz")
        assert len(results) == 0
    
    def test_transcript_to_gene_mapping(self, sample_transcriptome):
        """Test transcript to gene mapping functionality."""
        # Test getting gene from transcript
        gene_name = sample_transcriptome.get_gene_from_transcript("NM_001001234.1")
        assert gene_name == "Or9a"
        
        # Test with non-existent transcript
        gene_name = sample_transcriptome.get_gene_from_transcript("NonExistent")
        assert gene_name is None
    
    def test_get_transcript(self, sample_transcriptome):
        """Test transcript retrieval."""
        transcript = sample_transcriptome.get_transcript("NM_001001234.1")
        assert transcript is not None
        assert transcript.name == "NM_001001234.1"
        
        # Test with non-existent transcript
        transcript = sample_transcriptome.get_transcript("NonExistent")
        assert transcript is None


class TestGene:
    """Test Gene class functionality."""
    
    def test_initialization(self):
        """Test gene initialization."""
        gene = Gene("FBgn0000001", "Or9a")
        
        assert gene.id == "FBgn0000001"
        assert gene.name == "Or9a"
        assert isinstance(gene.transcripts, list)
        assert len(gene.transcripts) == 0
        assert gene.chromosome == ""
        assert gene.strand == ""
    
    def test_add_transcript(self, sample_transcript):
        """Test adding transcripts to gene."""
        gene = Gene("FBgn0000001", "Or9a")
        
        gene.add_transcript(sample_transcript)
        
        assert len(gene.transcripts) == 1
        assert gene.transcripts[0] == sample_transcript
    
    def test_get_transcript_longest_cds(self):
        """Test getting transcript with longest CDS."""
        gene = Gene("FBgn0000001", "TestGene")
        
        # Create transcripts with different CDS lengths
        transcript1 = Transcript("T1")
        transcript1.cds_length = 100
        
        transcript2 = Transcript("T2")
        transcript2.cds_length = 200  # Longest
        
        transcript3 = Transcript("T3")
        transcript3.cds_length = 150
        
        gene.add_transcript(transcript1)
        gene.add_transcript(transcript2)
        gene.add_transcript(transcript3)
        
        longest = gene.get_transcript_longest_cds()
        assert longest == transcript2
        assert longest.cds_length == 200
    
    def test_get_transcript_longest_bounds(self):
        """Test getting transcript with longest genomic bounds."""
        gene = Gene("FBgn0000001", "TestGene")
        
        # Create transcripts with different bounds
        transcript1 = Transcript("T1")
        transcript1.position = (1000, 1100)  # 100bp span
        
        transcript2 = Transcript("T2")
        transcript2.position = (2000, 2300)  # 300bp span (longest)
        
        transcript3 = Transcript("T3")
        transcript3.position = (3000, 3200)  # 200bp span
        
        # Add mock exons to make get_bounds() work
        for transcript in [transcript1, transcript2, transcript3]:
            exon = Exon(transcript.name, "ATCG", "chr1", transcript.position, "+")
            transcript.exons = [exon]
        
        gene.add_transcript(transcript1)
        gene.add_transcript(transcript2)
        gene.add_transcript(transcript3)
        
        longest = gene.get_transcript_longest_bounds()
        assert longest == transcript2
    
    def test_empty_gene_error(self):
        """Test error handling for genes without transcripts."""
        gene = Gene("FBgn0000001", "TestGene")
        
        with pytest.raises(ValueError, match="No transcripts found"):
            gene.get_transcript_longest_cds()
        
        with pytest.raises(ValueError, match="No transcripts found"):
            gene.get_transcript_longest_bounds()
    
    def test_string_representation(self):
        """Test gene string representation."""
        gene = Gene("FBgn0000001", "Or9a")
        gene.chromosome = "chr2L"
        gene.strand = "+"
        
        gene_str = str(gene)
        assert "Or9a" in gene_str
        assert "FBgn0000001" in gene_str
        assert "chr2L" in gene_str
        assert "+" in gene_str


class TestTranscript:
    """Test Transcript class functionality."""
    
    def test_initialization(self):
        """Test transcript initialization."""
        transcript = Transcript("NM_001001234.1")
        
        assert transcript.name == "NM_001001234.1"
        assert transcript.cds_sequence == ""
        assert transcript.cds_length == 0
        assert transcript.mrna_sequence == ""
        assert transcript.dna_sequence == ""
        assert transcript.chromosome == ""
        assert transcript.strand == ""
        assert transcript.position is None
        assert isinstance(transcript.utrs, list)
        assert isinstance(transcript.exons, list)
        assert isinstance(transcript.cds, list)
        assert isinstance(transcript.introns, list)
    
    def test_get_bounds(self, sample_transcript):
        """Test getting transcript bounds."""
        # Add some features to transcript
        exon1 = Exon("T1", "ATCG", "chr1", (1000, 1100), "+")
        exon2 = Exon("T1", "GCTA", "chr1", (1200, 1300), "+")
        utr = UTR("T1", "AAAA", "chr1", (900, 1000), "+", "5prime")
        
        sample_transcript.exons = [exon1, exon2]
        sample_transcript.utrs = [utr]
        
        bounds = sample_transcript.get_bounds()
        assert bounds is not None
        assert len(bounds) == 2
        assert bounds[0] == 900  # Start of UTR
        assert bounds[1] == 1300  # End of exon2
    
    def test_get_bounds_empty(self):
        """Test bounds with no features."""
        transcript = Transcript("T1")
        bounds = transcript.get_bounds()
        assert bounds is None
    
    def test_string_representation(self, sample_transcript):
        """Test transcript string representation."""
        transcript_str = str(sample_transcript)
        assert sample_transcript.name in transcript_str
        assert str(sample_transcript.cds_length) in transcript_str
        assert sample_transcript.chromosome in transcript_str


class TestGenomicFeatures:
    """Test genomic feature classes (Exon, CDS, UTR, Intron)."""
    
    def test_exon_creation(self):
        """Test Exon class."""
        exon = Exon("transcript1", "ATCGATCG", "chr2L", (1000, 1008), "+")
        
        assert exon.transcript == "transcript1"
        assert exon.sequence == "ATCGATCG"
        assert exon.chromosome == "chr2L"
        assert exon.position == (1000, 1008)
        assert exon.strand == "+"
        
        # Test string representation
        exon_str = str(exon)
        assert "Exon" in exon_str
        assert "transcript1" in exon_str
        assert "chr2L" in exon_str
        assert "1000-1008" in exon_str
    
    def test_cds_creation(self):
        """Test CDS class."""
        cds = CDS("transcript1", "ATGAAATAG", "chr2L", (1050, 1059), "+")
        
        assert cds.transcript == "transcript1"
        assert cds.sequence == "ATGAAATAG"
        assert cds.chromosome == "chr2L"
        assert cds.position == (1050, 1059)
        assert cds.strand == "+"
        
        # Test string representation
        cds_str = str(cds)
        assert "CDS" in cds_str
        assert "transcript1" in cds_str
        assert "chr2L" in cds_str
    
    def test_utr_creation(self):
        """Test UTR class."""
        utr = UTR("transcript1", "AAAAAAA", "chr2L", (900, 907), "+", "5prime")
        
        assert utr.transcript == "transcript1"
        assert utr.sequence == "AAAAAAA"
        assert utr.chromosome == "chr2L"
        assert utr.position == (900, 907)
        assert utr.strand == "+"
        assert utr.utr_type == "5prime"
        
        # Test string representation
        utr_str = str(utr)
        assert "UTR" in utr_str
        assert "5prime" in utr_str
        assert "transcript1" in utr_str
    
    def test_intron_creation(self):
        """Test Intron class."""
        intron = Intron("transcript1", "GTAG", "chr2L", (1100, 1104), "+")
        
        assert intron.transcript == "transcript1"
        assert intron.sequence == "GTAG"
        assert intron.chromosome == "chr2L"
        assert intron.position == (1100, 1104)
        assert intron.strand == "+"
        
        # Test string representation
        intron_str = str(intron)
        assert "Intron" in intron_str
        assert "transcript1" in intron_str


class TestIntegration:
    """Integration tests for transcriptomics module."""
    
    def test_complete_gene_structure(self):
        """Test creating a complete gene with all features."""
        # Create gene
        gene = Gene("FBgn0001", "TestGene")
        gene.chromosome = "chr2L"
        gene.strand = "+"
        
        # Create transcript
        transcript = Transcript("NM_001")
        transcript.chromosome = "chr2L"
        transcript.strand = "+"
        transcript.position = (1000, 2000)
        transcript.mrna_sequence = "A" * 500
        transcript.cds_sequence = "A" * 300
        transcript.cds_length = 300
        
        # Add features
        utr5 = UTR("NM_001", "A" * 50, "chr2L", (1000, 1050), "+", "5prime")
        exon1 = Exon("NM_001", "A" * 100, "chr2L", (1000, 1100), "+")
        intron1 = Intron("NM_001", "GT" + "A" * 96 + "AG", "chr2L", (1100, 1200), "+")
        exon2 = Exon("NM_001", "A" * 100, "chr2L", (1200, 1300), "+")
        cds1 = CDS("NM_001", "A" * 75, "chr2L", (1025, 1100), "+")
        cds2 = CDS("NM_001", "A" * 100, "chr2L", (1200, 1300), "+")
        utr3 = UTR("NM_001", "A" * 50, "chr2L", (1950, 2000), "+", "3prime")
        
        transcript.utrs = [utr5, utr3]
        transcript.exons = [exon1, exon2]
        transcript.introns = [intron1]
        transcript.cds = [cds1, cds2]
        
        gene.add_transcript(transcript)
        
        # Test structure
        assert len(gene.transcripts) == 1
        assert gene.get_transcript_longest_cds() == transcript
        assert gene.get_transcript_longest_bounds() == transcript
        
        # Test bounds calculation
        bounds = transcript.get_bounds()
        assert bounds[0] == 1000  # Start of first feature
        assert bounds[1] == 2000  # End of last feature
    
    def test_transcriptome_gene_lookup(self):
        """Test transcriptome gene and transcript lookup."""
        transcriptome = Transcriptome()
        
        # Create multiple genes with transcripts
        for i in range(3):
            gene = Gene(f"FBgn000{i}", f"Gene{i}")
            
            for j in range(2):  # 2 transcripts per gene
                transcript = Transcript(f"NM_00{i}00{j}")
                transcript.mrna_sequence = "ATCG" * 25
                gene.add_transcript(transcript)
            
            transcriptome.add_gene(gene)
        
        # Test gene retrieval
        gene0 = transcriptome.get_gene("Gene0")
        assert gene0 is not None
        assert len(gene0.transcripts) == 2
        
        # Test transcript lookup
        transcript = transcriptome.get_transcript("NM_001000")
        assert transcript is not None
        assert transcript.name == "NM_001000"
        
        # Test gene from transcript lookup
        gene_name = transcriptome.get_gene_from_transcript("NM_002001")
        assert gene_name == "Gene2"
        
        # Test search
        results = transcriptome.search_genes("Gene")
        assert len(results) == 3
    
    def test_transcript_to_gene_caching(self, sample_transcriptome):
        """Test that transcript-to-gene mapping is cached properly."""
        # First call should build cache
        gene_name1 = sample_transcriptome.get_gene_from_transcript("NM_001001234.1")
        assert len(sample_transcriptome.transcripts_to_genes) > 0
        
        # Second call should use cache
        gene_name2 = sample_transcriptome.get_gene_from_transcript("NM_001001234.1")
        assert gene_name1 == gene_name2
        
        # Adding new gene should clear cache
        new_gene = Gene("FBgn999", "NewGene")
        sample_transcriptome.add_gene(new_gene)
        assert len(sample_transcriptome.transcripts_to_genes) == 0