"""
Tests for hcrfish.hcr module.

This module tests the core HCR probe design functionality including
sequence manipulation, amplifier sequences, and probe generation.
"""

import pytest
import numpy as np
from typing import List, Tuple

from hcrfish.hcr.utils import (
    reverse_complement, 
    get_amplifier, 
    design_hcr_probes,
    reverse_string
)
from tests.conftest import VALID_AMPLIFIERS, INVALID_AMPLIFIERS, SAMPLE_SEQUENCES


class TestReverseComplement:
    """Test reverse complement function."""
    
    def test_basic_sequences(self):
        """Test basic DNA sequences."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GGGG") == "CCCC"
        assert reverse_complement("ATGC") == "GCAT"
    
    def test_case_sensitivity(self):
        """Test both upper and lower case."""
        assert reverse_complement("atcg") == "cgat"
        assert reverse_complement("ATcg") == "cgAT"
        assert reverse_complement("AtCg") == "cGaT"
    
    def test_special_characters(self):
        """Test handling of gaps and ambiguous bases."""
        assert reverse_complement("AT-CG") == "CG-AT"
        assert reverse_complement("ATNG") == "CNAT"
        assert reverse_complement("AT-N-GC") == "GC-N-AT"
    
    def test_empty_sequence(self):
        """Test empty sequence."""
        assert reverse_complement("") == ""
    
    def test_single_base(self):
        """Test single base sequences."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("G") == "C"
        assert reverse_complement("C") == "G"
    
    def test_palindromic_sequences(self):
        """Test palindromic sequences."""
        assert reverse_complement("ATAT") == "ATAT"
        assert reverse_complement("GCGC") == "GCGC"


class TestGetAmplifier:
    """Test amplifier sequence retrieval."""
    
    @pytest.mark.parametrize("amplifier", VALID_AMPLIFIERS)
    def test_valid_amplifiers(self, amplifier):
        """Test all valid amplifier types."""
        upspc, dnspc, up, dn = get_amplifier(amplifier)
        
        # Check that all components are strings
        assert isinstance(upspc, str)
        assert isinstance(dnspc, str) 
        assert isinstance(up, str)
        assert isinstance(dn, str)
        
        # Check spacer lengths (should be 2bp)
        assert len(upspc) == 2
        assert len(dnspc) == 2
        
        # Check initiator lengths (should be 18bp for HCR v3.0)
        assert len(up) == 18
        assert len(dn) == 18
        
        # Check that sequences contain only valid DNA bases (including lowercase)
        valid_bases = set("ATCGatcg")
        assert set(upspc).issubset(valid_bases)
        assert set(dnspc).issubset(valid_bases)
        assert set(up).issubset(valid_bases)
        assert set(dn).issubset(valid_bases)
    
    @pytest.mark.parametrize("invalid_amp", INVALID_AMPLIFIERS)
    def test_invalid_amplifiers(self, invalid_amp):
        """Test invalid amplifier types raise ValueError."""
        with pytest.raises(ValueError, match="Unsupported amplifier"):
            get_amplifier(invalid_amp)
    
    def test_amplifier_uniqueness(self):
        """Test that different amplifiers have different sequences."""
        amplifier_seqs = {}
        
        for amp in VALID_AMPLIFIERS:
            seqs = get_amplifier(amp)
            amplifier_seqs[amp] = seqs
        
        # Check that upstream initiators are unique
        upstream_seqs = [seqs[2] for seqs in amplifier_seqs.values()]
        assert len(set(upstream_seqs)) == len(VALID_AMPLIFIERS)
        
        # Check that downstream initiators are unique  
        downstream_seqs = [seqs[3] for seqs in amplifier_seqs.values()]
        assert len(set(downstream_seqs)) == len(VALID_AMPLIFIERS)
    
    def test_specific_amplifier_sequences(self):
        """Test specific known amplifier sequences."""
        # Test B3 (commonly used)
        upspc, dnspc, up, dn = get_amplifier("B3")
        assert upspc == "tt"
        assert dnspc == "tt"
        assert up == "GTCCCTGCCTCTATATCT"
        assert dn == "CCACTCAACTTTAACCCG"
        
        # Test B1
        upspc, dnspc, up, dn = get_amplifier("B1")
        assert upspc == "aa"
        assert dnspc == "ta"
        assert up == "GAGGAGGGCAGCAAACGG"
        assert dn == "GAAGAGTCTTCCTTTACG"


class TestDesignHcrProbes:
    """Test HCR probe design function."""
    
    def test_basic_probe_design(self, sample_mrna_sequence):
        """Test basic probe design functionality."""
        probes, regions, positions = design_hcr_probes(sample_mrna_sequence, "B3")
        
        # Should return lists
        assert isinstance(probes, list)
        assert isinstance(regions, list)
        assert isinstance(positions, list)
        
        # All lists should have same length
        assert len(probes) == len(regions) == len(positions)
        
        # Should generate at least one probe for a 200bp sequence
        assert len(probes) > 0
    
    def test_probe_structure(self, sample_mrna_sequence):
        """Test that probes have correct structure."""
        probes, regions, positions = design_hcr_probes(sample_mrna_sequence, "B3")
        
        if len(probes) > 0:
            probe_pair = probes[0]
            region = regions[0]
            position = positions[0]
            
            # Each probe pair should have exactly 2 probes
            assert len(probe_pair) == 2
            upstream_probe, downstream_probe = probe_pair
            
            # Get B3 amplifier sequences for comparison
            upspc, dnspc, up, dn = get_amplifier("B3")
            
            # Check probe structure
            assert upstream_probe.startswith(up)  # Should start with upstream initiator
            assert downstream_probe.endswith(dn)  # Should end with downstream initiator
            
            # Check region length (should be 52bp)
            assert len(region) == 52
            
            # Check position format
            assert len(position) == 2
            assert position[1] > position[0]  # End > start
            assert position[1] - position[0] == 52  # Correct span
    
    def test_gc_content_filtering(self):
        """Test that extreme GC content sequences are filtered."""
        # High GC sequence (should be filtered)
        high_gc_seq = "G" * 100 + "C" * 100  # 100% GC
        probes_hgc, _, _ = design_hcr_probes(high_gc_seq, "B3")
        
        # Low GC sequence (should be filtered)
        low_gc_seq = "A" * 100 + "T" * 100   # 0% GC
        probes_lgc, _, _ = design_hcr_probes(low_gc_seq, "B3")
        
        # Good GC sequence
        good_gc_seq = "ATCG" * 50  # 50% GC
        probes_good, _, _ = design_hcr_probes(good_gc_seq, "B3")
        
        # Good sequence should generate more probes than extreme sequences
        assert len(probes_good) >= len(probes_hgc)
        assert len(probes_good) >= len(probes_lgc)
    
    def test_homopolymer_filtering(self):
        """Test that sequences with long homopolymers are filtered."""
        # Sequence with long homopolymer runs
        homopoly_seq = "ATCG" * 20 + "GGGGG" + "ATCG" * 20  # 5 Gs in a row
        probes_homo, _, _ = design_hcr_probes(homopoly_seq, "B3")
        
        # Good sequence without homopolymers
        good_seq = "ATCG" * 50
        probes_good, _, _ = design_hcr_probes(good_seq, "B3")
        
        # Good sequence should generate more probes
        assert len(probes_good) >= len(probes_homo)
    
    def test_gap_filtering(self):
        """Test that sequences with gaps are filtered."""
        # Sequence with gaps (masked regions)
        gap_seq = "ATCG" * 20 + "-" * 10 + "ATCG" * 20
        probes_gap, _, _ = design_hcr_probes(gap_seq, "B3")
        
        # Good sequence without gaps
        good_seq = "ATCG" * 50
        probes_good, _, _ = design_hcr_probes(good_seq, "B3")
        
        # Good sequence should generate more probes
        assert len(probes_good) >= len(probes_gap)
    
    def test_short_sequence(self):
        """Test behavior with sequences too short for probe design."""
        short_seq = "ATCG"  # Only 4bp
        probes, regions, positions = design_hcr_probes(short_seq, "B3")
        
        # Should return empty lists
        assert len(probes) == 0
        assert len(regions) == 0
        assert len(positions) == 0
    
    def test_minimum_length_sequence(self):
        """Test with minimum length sequence (52bp)."""
        min_seq = "ATCG" * 13  # 52bp exactly
        probes, regions, positions = design_hcr_probes(min_seq, "B3")
        
        # Should generate exactly 1 probe
        assert len(probes) == 1
        assert len(regions) == 1
        assert len(positions) == 1
    
    @pytest.mark.parametrize("amplifier", VALID_AMPLIFIERS)
    def test_all_amplifiers(self, sample_mrna_sequence, amplifier):
        """Test probe design with all amplifier types."""
        probes, regions, positions = design_hcr_probes(sample_mrna_sequence, amplifier)
        
        # Should work with all amplifiers
        assert isinstance(probes, list)
        assert isinstance(regions, list)
        assert isinstance(positions, list)
        
        if len(probes) > 0:
            # Check that probes contain correct amplifier sequences
            upspc, dnspc, up, dn = get_amplifier(amplifier)
            
            upstream_probe, downstream_probe = probes[0]
            assert upstream_probe.startswith(up)
            assert downstream_probe.endswith(dn)
    
    def test_custom_parameters(self):
        """Test probe design with custom GC and homopolymer parameters."""
        seq = "ATCG" * 50
        
        # Strict parameters
        probes_strict, _, _ = design_hcr_probes(seq, "B3", gc_min=0.4, gc_max=0.6, max_homopolymer=3)
        
        # Relaxed parameters  
        probes_relaxed, _, _ = design_hcr_probes(seq, "B3", gc_min=0.1, gc_max=0.9, max_homopolymer=10)
        
        # Relaxed should generate same or more probes
        assert len(probes_relaxed) >= len(probes_strict)


class TestReverseString:
    """Test reverse string utility function."""
    
    def test_basic_reversal(self):
        """Test basic string reversal."""
        assert reverse_string("ATCG") == "GCTA"
        assert reverse_string("12345") == "54321"
        assert reverse_string("hello") == "olleh"
    
    def test_empty_string(self):
        """Test empty string."""
        assert reverse_string("") == ""
    
    def test_single_character(self):
        """Test single character."""
        assert reverse_string("A") == "A"
        assert reverse_string("1") == "1"
    
    def test_palindromes(self):
        """Test palindromic strings."""
        assert reverse_string("ATAT") == "TATA"
        assert reverse_string("12321") == "12321"
    
    def test_comparison_with_builtin(self):
        """Test that results match Python's built-in reversal."""
        test_strings = ["ATCG", "hello", "12345", "A", "", "ATAT"] 
        
        for s in test_strings:
            assert reverse_string(s) == s[::-1]


class TestIntegration:
    """Integration tests combining multiple functions."""
    
    def test_probe_design_pipeline(self, sample_transcriptome):
        """Test complete probe design pipeline."""
        # Get a gene from transcriptome
        gene = sample_transcriptome.get_gene("Or9a")
        assert gene is not None
        
        # Get transcript
        transcript = gene.get_transcript_longest_bounds()
        assert transcript is not None
        
        # Design probes
        probes, regions, positions = design_hcr_probes(transcript.mrna_sequence, "B3")
        
        # Should generate probes
        assert len(probes) > 0
        
        # Verify probe quality
        for probe_pair, region, position in zip(probes, regions, positions):
            # Check probe pair structure
            assert len(probe_pair) == 2
            
            # Check region
            assert len(region) == 52
            assert set(region).issubset(set("ATCG-"))
            
            # Check positions
            assert len(position) == 2
            assert 0 <= position[0] < position[1] <= len(transcript.mrna_sequence)
    
    def test_amplifier_probe_compatibility(self):
        """Test that probes work with their intended amplifiers."""
        seq = "ATCGATCGATCGATCG" * 10  # 160bp sequence
        
        for amplifier in VALID_AMPLIFIERS:
            probes, _, _ = design_hcr_probes(seq, amplifier)
            
            if len(probes) > 0:
                # Get amplifier sequences
                upspc, dnspc, up, dn = get_amplifier(amplifier)
                
                # Check each probe pair
                for upstream_probe, downstream_probe in probes:
                    # Verify structure
                    assert upstream_probe.startswith(up)
                    assert upstream_probe[len(up):len(up)+2] == upspc
                    assert downstream_probe.endswith(dn)
                    assert downstream_probe[-len(dn)-2:-len(dn)] == dnspc
                    
                    # Verify lengths (18bp initiator + 2bp spacer + 25bp target = 45bp)
                    assert len(upstream_probe) == 45
                    assert len(downstream_probe) == 45