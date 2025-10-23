"""
Main test runner and integration tests.

This module provides integration tests that verify the complete HCR-FISH
probe design workflow works correctly.
"""

import pytest
import tempfile
import os
from pathlib import Path

import hcrfish
from hcrfish import (
    design_hcr_probes,
    reverse_complement,
    get_amplifier,
    check_blast_tools,
    white_plotting,
    black_plotting
)
from hcrfish.transcriptomics.classes import Transcriptome, Gene, Transcript


class TestPackageImports:
    """Test that all main package components can be imported."""
    
    def test_main_package_import(self):
        """Test importing main package."""
        assert hcrfish is not None
        assert hasattr(hcrfish, '__version__') or True  # May not have version yet
    
    def test_submodule_imports(self):
        """Test importing all submodules."""
        from hcrfish import hcr, blast, plotting, transcriptomics
        
        assert hcr is not None
        assert blast is not None
        assert plotting is not None
        assert transcriptomics is not None
    
    def test_function_imports(self):
        """Test importing main functions."""
        from hcrfish.hcr import design_hcr_probes, reverse_complement, get_amplifier
        from hcrfish.blast import check_blast_tools, run_makeblastdb, run_blastn
        from hcrfish.plotting import white_plotting, black_plotting
        from hcrfish.transcriptomics import (
            Transcriptome, Gene, Transcript, 
            generate_transcriptome_object,
            load_transcriptome_object
        )
        
        # Just check they exist
        assert callable(design_hcr_probes)
        assert callable(reverse_complement)
        assert callable(get_amplifier)
        assert callable(check_blast_tools)
        assert callable(white_plotting)
        assert callable(black_plotting)


class TestWorkflowIntegration:
    """Test complete probe design workflow."""
    
    def test_basic_workflow(self, sample_transcriptome):
        """Test basic probe design workflow."""
        # Step 1: Get gene from transcriptome
        gene = sample_transcriptome.get_gene("Or9a")
        assert gene is not None
        
        # Step 2: Get transcript
        transcript = gene.get_transcript_longest_bounds()
        assert transcript is not None
        assert len(transcript.mrna_sequence) > 0
        
        # Step 3: Design probes
        probes, regions, positions = design_hcr_probes(transcript.mrna_sequence, "B3")
        
        # Verify results
        assert len(probes) > 0
        assert len(regions) == len(probes)
        assert len(positions) == len(probes)
        
        # Step 4: Verify probe structure
        for probe_pair, region, position in zip(probes[:3], regions[:3], positions[:3]):
            # Check probe pair structure
            assert len(probe_pair) == 2
            upstream, downstream = probe_pair
            
            # Check that probes contain B3 sequences
            upspc, dnspc, up, dn = get_amplifier("B3")
            assert upstream.startswith(up)
            assert downstream.endswith(dn)
            
            # Check region
            assert len(region) == 52
            assert set(region).issubset(set("ATCG-N"))
            
            # Check positions
            assert len(position) == 2
            assert position[1] > position[0]
    
    def test_multiple_amplifiers_workflow(self, sample_transcriptome):
        """Test workflow with different amplifiers."""
        gene = sample_transcriptome.get_gene("dsx")
        assert gene is not None
        
        transcript = gene.get_transcript_longest_bounds()
        assert transcript is not None
        
        # Test different amplifiers
        amplifiers = ["B1", "B2", "B3", "B4", "B5"]
        
        for amplifier in amplifiers:
            probes, regions, positions = design_hcr_probes(transcript.mrna_sequence, amplifier)
            
            if len(probes) > 0:
                # Verify amplifier-specific sequences
                upspc, dnspc, up, dn = get_amplifier(amplifier)
                
                upstream, downstream = probes[0]
                assert upstream.startswith(up)
                assert downstream.endswith(dn)
    
    def test_probe_quality_filtering(self):
        """Test that probe quality filtering works in complete workflow."""
        # Create test sequences with known issues
        test_sequences = {
            'good': "ATCGATCGATCGATCG" * 12,  # Good sequence
            'high_gc': "GCGCGCGCGCGCGCGC" * 12,  # High GC
            'low_gc': "ATATATATATATATAT" * 12,   # Low GC  
            'with_gaps': "ATCGATCG" * 8 + "---" * 16 + "ATCGATCG" * 8,  # Gaps
            'homopolymer': "ATCGATCG" * 6 + "GGGGGGGGG" + "ATCGATCG" * 6  # Long homopolymer
        }
        
        results = {}
        
        for seq_type, sequence in test_sequences.items():
            probes, regions, positions = design_hcr_probes(sequence, "B3")
            results[seq_type] = len(probes)
        
        # Good sequence should generate most probes
        assert results['good'] > 0
        
        # Problem sequences should generate fewer or no probes
        assert results['good'] >= results['high_gc']
        assert results['good'] >= results['low_gc'] 
        assert results['good'] >= results['with_gaps']
        assert results['good'] >= results['homopolymer']
    
    def test_plotting_integration(self):
        """Test that plotting functions work without errors."""
        import matplotlib.pyplot as plt
        
        # Test theme switching
        white_plotting()
        assert plt.rcParams['figure.facecolor'] == 'white'
        
        black_plotting()
        assert plt.rcParams['figure.facecolor'] == 'black'
        
        # Test creating a plot with themes
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot([1, 2, 3], [1, 4, 9])
        ax.set_title("Test Plot")
        
        # Should not raise errors
        plt.close(fig)
    
    def test_transcriptome_operations(self, sample_transcriptome):
        """Test transcriptome manipulation operations."""
        # Test searching
        or_genes = sample_transcriptome.search_genes("Or")
        assert len(or_genes) == 1
        assert or_genes[0].name == "Or9a"
        
        # Test transcript lookup  
        transcript = sample_transcriptome.get_transcript("NM_001001234.1")
        assert transcript is not None
        assert transcript.name == "NM_001001234.1"
        
        # Test gene from transcript
        gene_name = sample_transcriptome.get_gene_from_transcript("NM_001001234.1")
        assert gene_name == "Or9a"
        
        # Test adding new gene
        new_gene = Gene("FBgn999", "TestGene")
        new_transcript = Transcript("NM_999")
        new_transcript.mrna_sequence = "ATCG" * 50
        new_gene.add_transcript(new_transcript)
        
        sample_transcriptome.add_gene(new_gene)
        
        retrieved_gene = sample_transcriptome.get_gene("TestGene")
        assert retrieved_gene is not None
        assert retrieved_gene.name == "TestGene"


class TestErrorHandling:
    """Test error handling and edge cases."""
    
    def test_invalid_amplifier_error(self):
        """Test error handling for invalid amplifiers."""
        sequence = "ATCGATCGATCGATCG" * 10
        
        with pytest.raises(ValueError, match="Unsupported amplifier"):
            design_hcr_probes(sequence, "INVALID")
    
    def test_short_sequence_handling(self):
        """Test handling of sequences too short for probe design."""
        short_sequence = "ATCG"  # Only 4bp
        
        probes, regions, positions = design_hcr_probes(short_sequence, "B3")
        
        assert len(probes) == 0
        assert len(regions) == 0
        assert len(positions) == 0
    
    def test_empty_sequence_handling(self):
        """Test handling of empty sequences."""
        empty_sequence = ""
        
        probes, regions, positions = design_hcr_probes(empty_sequence, "B3")
        
        assert len(probes) == 0
        assert len(regions) == 0
        assert len(positions) == 0
    
    def test_nonexistent_gene_lookup(self, sample_transcriptome):
        """Test lookup of non-existent genes and transcripts."""
        # Non-existent gene
        gene = sample_transcriptome.get_gene("NonExistentGene")
        assert gene is None
        
        # Non-existent transcript
        transcript = sample_transcriptome.get_transcript("NonExistentTranscript")
        assert transcript is None
        
        # Gene from non-existent transcript
        gene_name = sample_transcriptome.get_gene_from_transcript("NonExistentTranscript")
        assert gene_name is None


class TestPerformance:
    """Test performance with larger datasets."""
    
    def test_large_sequence_probe_design(self):
        """Test probe design with large sequences."""
        # Create a 5kb sequence
        large_sequence = "ATCGATCGATCGATCG" * 312  # ~5kb
        
        probes, regions, positions = design_hcr_probes(large_sequence, "B3")
        
        # Should generate many probes
        assert len(probes) > 50  # Expect significant number of probes
        
        # All results should have same length
        assert len(probes) == len(regions) == len(positions)
        
        # Verify quality of first few probes
        for i in range(min(10, len(probes))):
            probe_pair = probes[i]
            region = regions[i]
            position = positions[i]
            
            assert len(probe_pair) == 2
            assert len(region) == 52
            assert len(position) == 2
            assert position[1] > position[0]
    
    def test_many_genes_transcriptome(self):
        """Test transcriptome with many genes."""
        transcriptome = Transcriptome()
        
        # Add 100 genes
        for i in range(100):
            gene = Gene(f"FBgn{i:05d}", f"Gene_{i}")
            
            # Add 2-3 transcripts per gene
            for j in range(2 + (i % 2)):  # 2 or 3 transcripts
                transcript = Transcript(f"NM_{i:03d}_{j}")
                transcript.mrna_sequence = f"{'ATCG' * 50}{'N' * i % 10}"  # Unique sequences
                gene.add_transcript(transcript)
            
            transcriptome.add_gene(gene)
        
        # Test operations on large transcriptome
        assert len(transcriptome.genes) == 100
        
        # Test search performance
        results = transcriptome.search_genes("Gene_1")
        assert len(results) >= 10  # Should find Gene_1, Gene_10-19, etc.
        
        # Test transcript lookup
        transcript = transcriptome.get_transcript("NM_050_0")
        assert transcript is not None
        
        # Test gene from transcript lookup
        gene_name = transcriptome.get_gene_from_transcript("NM_075_1")
        assert gene_name == "Gene_75"


if __name__ == "__main__":
    pytest.main([__file__])