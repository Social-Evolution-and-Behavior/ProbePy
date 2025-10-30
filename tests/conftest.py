"""
Test fixtures and utilities shared across test modules.

This module provides common test data, fixtures, and helper functions
used throughout the test suite.
"""

import pytest
import tempfile
import os
from pathlib import Path
from typing import Dict, Any

import numpy as np
from probepy.transcriptomics.classes import Transcriptome, Gene, Transcript, Exon, CDS, UTR


@pytest.fixture
def sample_dna_sequence():
    """Provide sample DNA sequence for testing."""
    return "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"


@pytest.fixture  
def sample_mrna_sequence():
    """Provide sample mRNA sequence suitable for probe design."""
    # 200bp sequence with good properties for probe design
    return ("ATGAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
            "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
            "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATA")


@pytest.fixture
def sample_transcript():
    """Create a sample transcript for testing."""
    transcript = Transcript("NM_001001234.1")
    transcript.mrna_sequence = ("ATGAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
                               "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                               "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
                               "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATA")
    transcript.cds_sequence = transcript.mrna_sequence[3:-3]  # Remove start/stop
    transcript.cds_length = len(transcript.cds_sequence)
    transcript.chromosome = "chr2L"
    transcript.strand = "+"
    transcript.position = (1000, 1200)
    return transcript


@pytest.fixture
def sample_gene():
    """Create a sample gene with transcript for testing."""
    gene = Gene("FBgn0000001", "Or9a")
    gene.chromosome = "chr2L"
    gene.strand = "+"
    
    # Add sample transcript
    transcript = Transcript("NM_001001234.1")
    transcript.mrna_sequence = ("ATGAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
                               "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                               "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
                               "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATA")
    transcript.cds_length = 194
    transcript.chromosome = "chr2L"
    transcript.strand = "+"
    transcript.position = (1000, 1200)
    
    # Add exons
    exon1 = Exon("NM_001001234.1", transcript.mrna_sequence[:100], "chr2L", (1000, 1100), "+")
    exon2 = Exon("NM_001001234.1", transcript.mrna_sequence[100:], "chr2L", (1100, 1200), "+")
    transcript.exons = [exon1, exon2]
    
    gene.add_transcript(transcript)
    return gene


@pytest.fixture
def sample_transcriptome():
    """Create a sample transcriptome with multiple genes."""
    transcriptome = Transcriptome()
    
    # Add first gene
    gene1 = Gene("FBgn0000001", "Or9a")
    gene1.chromosome = "chr2L"
    gene1.strand = "+"
    
    transcript1 = Transcript("NM_001001234.1")
    transcript1.mrna_sequence = ("ATGAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
                                "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                                "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
                                "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATA")
    transcript1.cds_length = 194
    transcript1.chromosome = "chr2L"
    transcript1.strand = "+"
    transcript1.position = (1000, 1200)
    gene1.add_transcript(transcript1)
    
    # Add second gene
    gene2 = Gene("FBgn0000002", "dsx")
    gene2.chromosome = "chr3R"
    gene2.strand = "-"
    
    transcript2 = Transcript("NM_001001235.1")
    transcript2.mrna_sequence = ("GCGCGCATGAAATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
                                "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
                                "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
                                "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGCGCGC")
    transcript2.cds_length = 194
    transcript2.chromosome = "chr3R"
    transcript2.strand = "-"
    transcript2.position = (2000, 2200)
    gene2.add_transcript(transcript2)
    
    transcriptome.add_gene(gene1)
    transcriptome.add_gene(gene2)
    
    return transcriptome


@pytest.fixture
def temp_fasta_file():
    """Create a temporary FASTA file for testing."""
    content = """>test_sequence_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>test_sequence_2  
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(content)
        temp_path = f.name
    
    yield temp_path
    
    # Cleanup
    try:
        os.unlink(temp_path)
    except (OSError, FileNotFoundError):
        pass


@pytest.fixture
def temp_directory():
    """Create a temporary directory for testing."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    
    # Cleanup
    import shutil
    try:
        shutil.rmtree(temp_dir)
    except (OSError, FileNotFoundError):
        pass


@pytest.fixture
def simple_genome_file():
    """Create a simple genome FASTA file for testing."""
    content = """>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(content)
        temp_path = f.name
    
    yield temp_path
    
    # Cleanup
    try:
        os.unlink(temp_path)
    except (OSError, FileNotFoundError):
        pass


# Test data constants
VALID_AMPLIFIERS = ["B1", "B2", "B3", "B4", "B5"]
INVALID_AMPLIFIERS = ["B6", "A1", "invalid", ""]

SAMPLE_SEQUENCES = {
    'good_gc': "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    'high_gc': "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC", 
    'low_gc': "ATATATATATATATATATATATATATATATATATATATATATATAT",
    'with_gaps': "ATCGATCG---ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    'homopolymer': "ATCGATCGGGGGGATCGATCGATCGATCGATCGATCGATCGATCGAT"
}


class TestDataGenerator:
    """Helper class for generating test data."""
    
    @staticmethod
    def create_gene(gene_id: str, gene_name: str, 
                   transcript_count: int = 1) -> Gene:
        """Create a gene with specified number of transcripts."""
        gene = Gene(gene_id, gene_name)
        gene.chromosome = "chr2L"
        gene.strand = "+"
        
        for i in range(transcript_count):
            transcript = Transcript(f"NM_{gene_id}_{i+1}")
            transcript.mrna_sequence = SAMPLE_SEQUENCES['good_gc'] * 4  # 188bp
            transcript.cds_length = 180
            transcript.chromosome = gene.chromosome
            transcript.strand = gene.strand
            transcript.position = (1000 + i*100, 1200 + i*100)
            gene.add_transcript(transcript)
            
        return gene
    
    @staticmethod
    def create_random_sequence(length: int, gc_content: float = 0.5) -> str:
        """Generate random DNA sequence with specified GC content."""
        np.random.seed(42)  # For reproducible tests
        
        gc_bases = int(length * gc_content)
        at_bases = length - gc_bases
        
        bases = ['G', 'C'] * (gc_bases // 2) + ['A', 'T'] * (at_bases // 2)
        
        # Handle odd numbers
        if gc_bases % 2:
            bases.append('G')
        if at_bases % 2:
            bases.append('A')
            
        np.random.shuffle(bases)
        return ''.join(bases)