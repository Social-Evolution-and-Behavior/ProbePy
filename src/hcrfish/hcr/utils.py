"""
HCR-FISH probe design utilities.

This module contains core functions for designing HCR v3.0 compatible probes,
including sequence manipulation and amplifier sequence retrieval.
"""

from typing import Tuple, List, Optional
import numpy as np


def reverse_complement(sequence: str) -> str:
    """
    Get the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): Input DNA sequence (A, T, G, C, -)
        
    Returns:
        str: Reverse complement of the input sequence
        
    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("ATCG-N")
        'N-CGAT'
    """
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        '-': '-', 'N': 'N', 'n': 'n'
    }
    
    complement = ''
    for base in sequence:
        complement = complement_map.get(base, base) + complement
    return complement


def get_amplifier(amplifier: str) -> Tuple[str, str, str, str]:
    """
    Get HCR v3.0 amplifier sequences and spacers.
    
    Args:
        amplifier (str): Amplifier name (B1, B2, B3, B4, or B5)
        
    Returns:
        Tuple[str, str, str, str]: (upstream_spacer, downstream_spacer, 
                                   upstream_initiator, downstream_initiator)
        
    Raises:
        ValueError: If amplifier is not supported
        
    Examples:
        >>> upspc, dnspc, up, dn = get_amplifier("B3")
        >>> print(f"B3 upstream: {up}")
        B3 upstream: GTCCCTGCCTCTATATCT
    """
    amplifiers = {
        "B1": {
            "upstream_spacer": "aa",
            "downstream_spacer": "ta", 
            "upstream_initiator": "GAGGAGGGCAGCAAACGG",
            "downstream_initiator": "GAAGAGTCTTCCTTTACG"
        },
        "B2": {
            "upstream_spacer": "aa",
            "downstream_spacer": "aa",
            "upstream_initiator": "CCTCGTAAATCCTCATCA", 
            "downstream_initiator": "ATCATCCAGTAAACCGCC"
        },
        "B3": {
            "upstream_spacer": "tt",
            "downstream_spacer": "tt",
            "upstream_initiator": "GTCCCTGCCTCTATATCT",
            "downstream_initiator": "CCACTCAACTTTAACCCG"
        },
        "B4": {
            "upstream_spacer": "aa",
            "downstream_spacer": "at",
            "upstream_initiator": "CCTCAACCTACCTCCAAC",
            "downstream_initiator": "TCTCACCATATTCGCTTC"
        },
        "B5": {
            "upstream_spacer": "aa", 
            "downstream_spacer": "aa",
            "upstream_initiator": "CTCACTCCCAATCTCTAT",
            "downstream_initiator": "CTACCCTACAAATCCAAT"
        }
    }
    
    if amplifier not in amplifiers:
        raise ValueError(f"Unsupported amplifier: {amplifier}. "
                        f"Supported amplifiers: {list(amplifiers.keys())}")
    
    amp = amplifiers[amplifier]
    return (
        amp["upstream_spacer"],
        amp["downstream_spacer"], 
        amp["upstream_initiator"],
        amp["downstream_initiator"]
    )


def design_hcr_probes(sequence: str, amplifier: str, 
                      gc_min: float = 0.25, gc_max: float = 0.75,
                      max_homopolymer: int = 4) -> Tuple[List[List[str]], List[str], List[List[int]]]:
    """
    Design HCR-FISH probe pairs for a given target sequence.
    
    This function designs probe pairs compatible with HCR v3.0 amplifiers by:
    1. Working 3' to 5' along the reverse complement of the target
    2. Checking each potential 52bp probe region for quality issues
    3. Generating complementary 25nt probe pairs with appropriate spacers
    
    Args:
        sequence (str): Target mRNA sequence (mature transcript, no introns)
        amplifier (str): HCR amplifier type (B1, B2, B3, B4, or B5)
        gc_min (float, optional): Minimum GC content for probe segments. Defaults to 0.25.
        gc_max (float, optional): Maximum GC content for probe segments. Defaults to 0.75.
        max_homopolymer (int, optional): Maximum homopolymer run length. Defaults to 4.
        
    Returns:
        Tuple containing:
            - probe_pairs (List[List[str]]): List of [upstream_probe, downstream_probe] pairs
            - probe_regions (List[str]): List of 52bp target regions where probes bind
            - probe_positions (List[List[int]]): List of [start, end] positions in original sequence
            
    Examples:
        >>> sequence = "ATCGATCGATCG" * 10  # Example target sequence
        >>> probes, regions, positions = design_hcr_probes(sequence, "B3")
        >>> print(f"Designed {len(probes)} probe pairs")
        
    Notes:
        - Probes are 25 nucleotides each with a 2bp overlap region
        - Quality filters reject regions with gaps (-), extreme GC content, or long homopolymers
        - Probe spacing is optimized for HCR v3.0 performance (54bp between probe centers)
    """
    # Convert to reverse complement for 3' to 5' design
    sequence_rc = reverse_complement(sequence)
    
    # Get amplifier-specific sequences
    upstream_spacer, downstream_spacer, upstream_init, downstream_init = get_amplifier(amplifier)
    
    # Initialize output lists
    probe_pairs: List[List[str]] = []
    probe_regions: List[str] = []
    probe_positions: List[List[int]] = []
    
    # Start from 3' end and work towards 5'
    position = len(sequence_rc)
    probe_length = 52  # Total probe pair binding region
    
    while position >= probe_length:
        # Extract potential probe binding regions (25bp each)
        upstream_region = str(sequence_rc[position-25:position])
        downstream_region = str(sequence_rc[position-52:position-27])
        
        # Calculate GC content for each region
        gc_upstream = (upstream_region.count("G") + upstream_region.count("C")) / 25
        gc_downstream = (downstream_region.count("G") + downstream_region.count("C")) / 25
        
        # Quality control checks
        skip_position = False
        skip_amount = 1
        
        # Check for gaps (masked regions)
        if '-' in upstream_region or '-' in downstream_region:
            skip_position = True
            
        # Check GC content limits
        elif (gc_upstream > gc_max or gc_upstream < gc_min or 
              gc_downstream > gc_max or gc_downstream < gc_min):
            skip_position = True
            
        # Check for problematic homopolymers
        else:
            homopolymers = ['G' * (max_homopolymer + 1), 'C' * (max_homopolymer + 1),
                          'A' * (max_homopolymer + 1), 'T' * (max_homopolymer + 1)]
            
            for homopoly in homopolymers:
                if homopoly in upstream_region or homopoly in downstream_region:
                    skip_position = True
                    skip_amount = max_homopolymer + 1
                    break
        
        if skip_position:
            position -= skip_amount
            continue
            
        # Generate probe sequences with amplifier-specific sequences
        upstream_probe = upstream_init + upstream_spacer + upstream_region
        downstream_probe = downstream_region + downstream_spacer + downstream_init
        
        # Store results
        probe_pairs.append([upstream_probe, downstream_probe])
        
        # Store the target region (reverse complement back to original orientation)
        target_region = reverse_complement(sequence_rc[position-52:position])
        probe_regions.append(target_region)
        
        # Calculate positions in original sequence coordinates
        start_pos = len(sequence) - position
        end_pos = start_pos + probe_length
        probe_positions.append([start_pos, end_pos])
        
        # Move to next position (54bp spacing for optimal HCR performance)
        position -= 54
    
    return probe_pairs, probe_regions, probe_positions


def reverse_string(string: str) -> str:
    """
    Reverse a string using numpy indexing.
    
    Args:
        string (str): Input string to reverse
        
    Returns:
        str: Reversed string
        
    Examples:
        >>> reverse_string("ATCG")
        'GCTA'
        
    Note:
        This function is kept for backwards compatibility but str[::-1] is more efficient.
    """
    if not string:
        return string
    indices = np.linspace(len(string)-1, 0, len(string)).astype(int)
    return ''.join([string[i] for i in indices])