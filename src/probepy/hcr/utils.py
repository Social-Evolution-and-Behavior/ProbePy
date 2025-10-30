"""
HCR-FISH probe design utilities.

This module contains core functions for designing HCR v3.0 compatible probes,
including sequence manipulation and amplifier sequence retrieval.
"""

import re
import subprocess
from typing import Tuple, List, Optional, Any
import numpy as np
import os
import pandas as pd
from pygenomeviz import GenomeViz
import Bio.SeqIO as SeqIO
import matplotlib.pyplot as plt

from probepy.blast.install import check_blast_tools
from probepy.transcriptomics.classes import Gene, Transcriptome


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



def blast_gene(
    gene_name: str,
    transcriptome: Transcriptome,
    species_identifier: str,
    permitted_off_targets: Optional[List[str]] = None,
    length_thresh: int = 50,
    base_dir: str = "",
) -> None:
    """
    BLAST gene sequences against transcriptome databases to identify potential off-targets.
    
    This function performs BLAST searches against both mature mRNA (no introns) and 
    pre-mRNA (with introns) databases to identify regions of similarity that could 
    interfere with probe specificity.
    
    Args:
        gene_name (str): Name of the target gene
        transcriptome (Transcriptome): Transcriptome object with gene annotation
        species_identifier (str): Species identifier for database paths
        permitted_off_targets (Optional[List[str]], optional): List of permitted off-target keywords. Defaults to None.
        length_thresh (int, optional): Minimum alignment length threshold. Defaults to 50.
        base_dir (str, optional): Base directory for input/output files. Defaults to "".
        
    Raises:
        Exception: If BLAST+ tools are not available
        ValueError: If gene is None or lacks target_sequence attribute
        subprocess.CalledProcessError: If BLAST commands fail
        
    Notes:
        - Creates FASTA files in output/{species}/gene_seq_blast_input/
        - Outputs results to output/{species}/gene_seq_blast_output/
        - Uses ungapped BLAST with strict parameters for sensitivity
        - Searches both mature and pre-mRNA transcriptomes
    """

    # Validate BLAST tools availability
    check = check_blast_tools()
    if not check['blastn']['available']:
        raise Exception(
            "blastn not found. Please install BLAST+ tools and ensure they are in your PATH."
        )
    
    gene = transcriptome.get_gene(gene_name)

    # Confirm gene has a target sequence
    if not gene or not gene.target_sequence:
        raise ValueError("Gene not found or does not have a target_sequence.")
    
    if not hasattr(gene, 'target_sequence') or gene.target_sequence is None:
        raise ValueError("Gene does not have a target_sequence attribute.")
    
    # Define BLAST database paths
    transcriptome_no_introns_db = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_no_introns", "mRNA_no_introns")
    transcriptome_yes_introns_db = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_yes_introns", "mRNA_yes_introns")
    

    # Prepare input directory and FASTA file
    gene_seq_blast_input_dir = os.path.join(base_dir, "output", species_identifier, 'gene_seq_blast_input')
    os.makedirs(gene_seq_blast_input_dir, exist_ok=True)

    # Clear existing files in input directory
    for file in os.listdir(gene_seq_blast_input_dir):
        os.remove(os.path.join(gene_seq_blast_input_dir, file))

    # Write gene sequence to FASTA file
    input_fasta_path = os.path.join(gene_seq_blast_input_dir, f"{gene.name}.fasta")
    with open(input_fasta_path, 'w') as f:
        f.write(f">{gene.name}\n{gene.target_sequence}")

    # Prepare output directory
    gene_seq_blast_output_dir = os.path.join(base_dir, "output", species_identifier, 'gene_seq_blast_output')
    os.makedirs(gene_seq_blast_output_dir, exist_ok=True)

    # Clear existing files in output directory
    for file in os.listdir(gene_seq_blast_output_dir):
        os.remove(os.path.join(gene_seq_blast_output_dir, file))

    # BLAST parameters for sensitive detection
    blast_params = [
        "-task blastn",
        f"-query {input_fasta_path}",
        "-ungapped",
        "-word_size 15",
        "-strand plus",
        "-reward 1",
        "-penalty -5",
        "-dust no",
        "-soft_masking false",
        "-max_target_seqs 10000",
        "-outfmt '10 qseqid sseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore'",
        "-num_threads 4"
    ]

    # BLAST against mature mRNA database (no introns)
    output_path_no_introns = os.path.join(
        gene_seq_blast_output_dir, f"{gene.name}_blasted_no_introns.csv"
    )
    command_no_introns = " ".join([
        "blastn"
    ] + blast_params[1:] + [
        f"-db {transcriptome_no_introns_db}",
        f"-out {output_path_no_introns}"
    ])
    
    try:
        subprocess.run(command_no_introns, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"BLAST against no-introns database failed: {e}")

    # BLAST against pre-mRNA database (with introns)
    output_path_yes_introns = os.path.join(
        gene_seq_blast_output_dir, f"{gene.name}_blasted_yes_introns.csv"
    )
    command_yes_introns = " ".join([
        "blastn"
    ] + blast_params[1:] + [
        f"-db {transcriptome_yes_introns_db}",
        f"-out {output_path_yes_introns}"
    ])
    
    try:
        subprocess.run(command_yes_introns, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"BLAST against introns database failed: {e}")

    print(f"BLAST searches completed for {gene.name}. Results saved to {gene_seq_blast_output_dir}")
    
    if permitted_off_targets is None:
        permitted_off_targets = []
    
    # Define input and output directories
    gene_seq_blast_output_dir = os.path.join(base_dir, "output", species_identifier, 'gene_seq_blast_output')
    gene_seq_unique_regions_dir = os.path.join(base_dir, "output", species_identifier, 'gene_seq_unique_regions')
    
    # Create output directory
    os.makedirs(gene_seq_unique_regions_dir, exist_ok=True)

    # Clear existing files in output directory
    for file in os.listdir(gene_seq_unique_regions_dir):
        os.remove(os.path.join(gene_seq_unique_regions_dir, file))

    # Define column names for BLAST output format
    blast_columns = [
        'query_id', 'subject_id', 'subject_acc', 'percent_identity', 'length',
        'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end',
        'evalue', 'bitscore'
    ]

    # Load BLAST results from both databases
    try:
        # No introns results
        no_introns_path = os.path.join(gene_seq_blast_output_dir, f"{gene.name}_blasted_no_introns.csv")
        blast_results_no_introns = pd.read_csv(
            no_introns_path, header=None, names=blast_columns
        )
        blast_results_no_introns['source'] = 'no_introns'

        # With introns results  
        yes_introns_path = os.path.join(gene_seq_blast_output_dir, f"{gene.name}_blasted_yes_introns.csv")
        blast_results_yes_introns = pd.read_csv(
            yes_introns_path, header=None, names=blast_columns
        )
        blast_results_yes_introns['source'] = 'yes_introns'
        
    except FileNotFoundError as e:
        raise FileNotFoundError(f"BLAST result files not found: {e}")

    # Combine results from both databases
    blast_results = pd.concat([blast_results_no_introns, blast_results_yes_introns], axis=0)

    # Filter by alignment length threshold
    blast_results = blast_results[blast_results['length'] >= length_thresh].reset_index(drop=True)

    # Annotate with subject gene information
    blast_results['subject_gene_id'] = blast_results.apply(
        lambda x: (
            transcriptome.get_gene_from_transcript(x['subject_acc']) or
            transcriptome.get_gene_from_transcript(x['subject_id'])
        ), 
        axis=1
    )

    # Identify same-gene hits (self-matches)
    blast_results['same_gene'] = (
        blast_results['query_id'] == blast_results['subject_gene_id']
    )

    # Remove empty strings from permitted_off_targets
    permitted_off_targets = [term for term in permitted_off_targets if term]

    # Identify permitted off-targets
    blast_results['permitted_off_target'] = blast_results['subject_gene_id'].apply(
        lambda subject_id: any(keyword in subject_id for keyword in permitted_off_targets)
    )

    # Get target sequence for masking
    if not hasattr(gene, 'target_sequence') or gene.target_sequence is None:
        raise ValueError("Gene does not have a target_sequence attribute.")
    
    sequence = gene.target_sequence

    # Identify problematic off-target regions
    off_targets = blast_results[
        (blast_results['length'] >= length_thresh) &
        (~blast_results['same_gene']) &
        (~blast_results['permitted_off_target'])
    ]

    # Mask off-target regions with '-' characters
    for _, row in off_targets.iterrows():
        start_pos = int(row['q_start']) - 1  # Convert to 0-based indexing
        end_pos = int(row['q_end'])
        mask_length = end_pos - start_pos
        sequence = sequence[:start_pos] + '-' * mask_length + sequence[end_pos:]

    # Store results in gene object
    gene.unique_sequence = sequence
    gene.blast_results = blast_results
    gene.off_targets = off_targets

    # Export masked sequence to file
    output_fasta_path = os.path.join(gene_seq_unique_regions_dir, f"{gene.name}_unique.fasta")
    with open(output_fasta_path, 'w') as f:
        f.write(f">{gene.name}\n{sequence}")

    print(f"Unique regions have been annotated and exported to {gene_seq_unique_regions_dir}")
    print(f"Masked {len(off_targets)} off-target regions from {gene.name}")


def get_probes_IDT(
    gene_name: str,
    transcriptome: Transcriptome,
    species_identifier: str,
    amplifier: str = "B1",
    n_probes: int = 30,
    base_dir: str = "",
) -> None:
    """
    Design HCR-FISH probes and export them in IDT-compatible format.
    
    This function designs HCR-FISH probes for a gene's unique sequence regions
    and exports them as Excel files suitable for ordering from IDT (Integrated 
    DNA Technologies).
    
    Args:
        gene_name (str): Name of the target gene
        transcriptome (Transcriptome): Transcriptome object with gene annotation
        species_identifier (str): Species identifier for file organization
        amplifier (str, optional): HCR amplifier type (B1-B5). Defaults to "B1".
        n_probes (int, optional): Maximum number of probe pairs to select. Defaults to 30.
        base_dir (str, optional): Base directory for output files. Defaults to "".
            
    Raises:
        ValueError: If gene lacks unique_sequence attribute
        ValueError: If amplifier is not supported
        
    Notes:
        - Randomly selects n_probes if more are available
        - Exports probe sequences to output/{species}/IDT_sheets/
        - Exports probe binding regions to output/{species}/probe_binding_regions_sheets/
        - Uses current date in filenames
        - Stores probe data in gene.probes and gene.regions attributes
    """

    gene = transcriptome.get_gene(gene_name)

    if gene is None:
        raise ValueError(f"Gene '{gene_name}' not found in transcriptome.")
    
    # Validate gene has processed unique sequence
    if not hasattr(gene, 'unique_sequence') or gene.unique_sequence is None:
        raise ValueError(
            "Gene does not have a unique_sequence attribute. "
            "Please run process_gene_blast_results first."
        )

    # Design probes using the unique sequence
    try:
        probes, regions, positions = design_hcr_probes(gene.unique_sequence, amplifier)
    except ValueError as e:
        raise ValueError(f"Probe design failed: {e}")
    
    print(f"{len(probes)} probes designed for {gene.name} using amplifier {amplifier}")

    # Select subset of probes if needed
    np.random.seed(1)  # For reproducible selection
    if len(probes) <= n_probes:
        selected_indices = list(range(len(probes)))
        print(f"Using all {len(probes)} available probes for {gene.name}")
    else:
        selected_indices = np.random.choice(
            range(len(probes)), n_probes, replace=False
        ).tolist()
        print(f"Randomly selected {n_probes} probes from {len(probes)} available for {gene.name}")

    # Extract selected probes and regions
    selected_probes = [probes[i] for i in selected_indices]
    selected_regions = [regions[i] for i in selected_indices]

    # Store results in gene object
    gene.probes = selected_probes
    gene.regions = selected_regions

    # Prepare probe sequences for IDT format (flatten probe pairs)
    probe_sequences = []
    for probe_pair in selected_probes:
        probe_sequences.extend(probe_pair)

    # Generate timestamp for filenames
    timestamp = pd.Timestamp.now().strftime('%Y-%m-%d')

    # Export probe sequences for IDT ordering
    idt_output_dir = os.path.join(base_dir, "output", species_identifier, 'IDT_sheets')
    os.makedirs(idt_output_dir, exist_ok=True)
    
    idt_output_path = os.path.join(
        idt_output_dir, f"{gene.name}-{amplifier}-{timestamp}.xlsx"
    )
    
    # Create IDT-formatted DataFrame
    idt_df = pd.DataFrame({
        'Pool name': [f'{gene.name}-{amplifier}'] * len(probe_sequences),
        'Sequence': probe_sequences
    })
    
    try:
        idt_df.to_excel(idt_output_path, index=False)
        print(f"Exported {len(probe_sequences)} probe sequences to {idt_output_path}")
    except Exception as e:
        raise Exception(f"Failed to export IDT file: {e}")

    # Export probe binding regions for reference
    regions_output_dir = os.path.join(
        base_dir, "output", species_identifier, 'probe_binding_regions_sheets'
    )
    os.makedirs(regions_output_dir, exist_ok=True)
    
    regions_output_path = os.path.join(
        regions_output_dir, f"{gene.name}-{amplifier}-regions-{timestamp}.xlsx"
    )
    
    # Create probe regions DataFrame
    regions_df = pd.DataFrame({
        'Gene': [gene.name] * len(selected_probes),
        'Region': selected_regions,
        'Probe 1': [probe_pair[0] for probe_pair in selected_probes],
        'Probe 2': [probe_pair[1] for probe_pair in selected_probes]
    })
    
    try:
        regions_df.to_excel(regions_output_path, index=False)
        print(f"Exported probe binding regions to {regions_output_path}")
    except Exception as e:
        raise Exception(f"Failed to export regions file: {e}")


def get_probe_binding_regions_plot(
    gene_name: str,
    transcriptome: Transcriptome,
    base_dir: str,
    species_identifier: str,
    save: bool = True
) -> plt.Figure:
    """
    Generate and return a genomic visualization of probe binding regions.
    
    This function creates a publication-ready plot showing probe binding sites
    mapped to the gene's genomic structure, including exons, UTRs, and introns.
    
    Args:
        gene_name (str): Name of the target gene
        transcriptome (Transcriptome): Transcriptome object with gene annotation
        base_dir (str): Base directory for input/output files
        species_identifier (str): Species identifier for file organization
        save (bool, optional): Whether to save the plot to disk. Defaults to True.
        
    Returns:
        plt.Figure: Matplotlib figure object containing the genomic visualization
        
    Raises:
        ValueError: If gene lacks regions attribute
        FileNotFoundError: If genome FASTA file is not found
        Exception: If genomic visualization fails
        
    Notes:
        - Requires genome FASTA file at input/{species}/genome/genome_fasta.fa
        - Maps probe regions to both forward and reverse genomic strands
        - If save=True, exports high-resolution PNG to output/{species}/probe_regions_plot/
        - Returns matplotlib figure for display in notebooks
    """

    gene = transcriptome.get_gene(gene_name)
    if gene is None:
        raise ValueError(f"Gene '{gene_name}' not found in transcriptome.")
    
    # Validate gene has probe regions
    if not hasattr(gene, 'regions') or gene.regions is None:
        raise ValueError(
            "Gene does not have a regions attribute. "
            "Please run get_probes_IDT first."
        )
    
    regions = gene.regions

    # Load genome sequence
    genome_fasta_dir = os.path.join(base_dir, "input", species_identifier, "genome")

    # Look for a file ending with .fa or .fasta
    genome_fasta_path = None
    for file in os.listdir(genome_fasta_dir):
        if file.endswith('.fa') or file.endswith('.fasta'):
            genome_fasta_path = os.path.join(genome_fasta_dir, file)
            break

    if genome_fasta_path is None:
        raise FileNotFoundError(f"No genome FASTA file found in {genome_fasta_dir}")

     # Load genome sequences into a dictionary
    try:
        genome_seq = SeqIO.to_dict(SeqIO.parse(genome_fasta_path, "fasta"))
    except FileNotFoundError:
        raise FileNotFoundError(f"Genome FASTA file not found: {genome_fasta_path}")
    except Exception as e:
        raise Exception(f"Failed to load genome sequence: {e}")

    try:
        # Initialize genomic visualization
        gv = GenomeViz(track_align_type="center")
        gv.set_scale_bar(ymargin=0.5)

        # Get gene coordinates and chromosome
        transcript = gene.get_transcript_longest_bounds()
        chromosome = gene.chromosome
        bounds = transcript.get_bounds()
        min_start = int(bounds[0])
        max_end = int(bounds[1])
        
        # Create main genomic track
        track = gv.add_feature_track(chromosome, segments=(min_start, max_end))
        track.add_sublabel()

        # Add exons with directional arrows
        exon_bounds = [
            (int(exon.position[0]), int(exon.position[1])) 
            for exon in transcript.exons
        ]
        exon_bounds = sorted(exon_bounds, key=lambda x: x[0])
        strand = 1 if transcript.strand == '+' else -1
        
        track.add_exon_feature(
            exon_bounds, 
            strand, 
            plotstyle='arrow', 
            arrow_shaft_ratio=0.5, 
            patch_kws=dict(fc="green", ec="none", alpha=0.5), 
            intron_patch_kws=dict(ec="black", lw=2), 
            label=''
        )

        # Add UTR regions
        for utr in transcript.utrs:
            track.add_feature(
                int(utr.position[0]), 
                int(utr.position[1]), 
                strand, 
                plotstyle="box", 
                lw=0, 
                fc='grey', 
                alpha=0.5
            )

        # Extract genomic sequence for probe mapping
        if chromosome not in genome_seq:
            raise ValueError(f"Chromosome {chromosome} not found in genome sequence")
        
        forward_seq = str(genome_seq[chromosome].seq[min_start:max_end]).upper()
        reverse_seq = reverse_complement(forward_seq)

        # Map probe regions to forward strand
        for region in regions:
            forward_positions = [m.start() for m in re.finditer(region, forward_seq)]
            for position in forward_positions:
                track.add_feature(
                    min_start + position, 
                    min_start + position + len(region), 
                    1, 
                    plotstyle="box", 
                    label='', 
                    ec="none", 
                    fc="magenta", 
                    alpha=0.8
                )

        # Map probe regions to reverse strand
        for region in regions:
            reverse_positions = [m.start() for m in re.finditer(region, reverse_seq)]
            for position in reverse_positions:
                track.add_feature(
                    max_end - position - len(region), 
                    max_end - position, 
                    -1, 
                    plotstyle="box", 
                    label='', 
                    ec="none", 
                    fc="blue", 
                    alpha=0.8
                )

        # Generate plot
        fig = gv.plotfig()
        plt.title(f"{gene.name} HCR-FISH Probe Binding Sites", y=1.8, fontsize=16)
        
        # Save to disk if requested
        if save:
            output_dir = os.path.join(base_dir, "output", species_identifier, 'probe_regions_plot')
            os.makedirs(output_dir, exist_ok=True)
            
            output_path = os.path.join(output_dir, f"{gene.name}-probes.png")
            fig.savefig(output_path, bbox_inches='tight', dpi=300)
            
            print(f"Probe binding regions plot exported to {output_path}")
        
        return fig
        
    except Exception as e:
        plt.close('all')
        raise Exception(f"Failed to generate genomic visualization: {e}")

def check_probe_availability(
    gene_name: str,
    transcriptome: Transcriptome,
    species_identifier: str,
    base_dir: str = "",
    permitted_off_targets: Optional[List[str]] = None,
) -> int:
    """
    Check the number of HCR-FISH probes available for a given gene.
    
    This function performs the complete probe design pipeline to determine
    how many high-quality probe pairs can be generated for a target gene.
    It includes BLAST analysis to identify unique regions and probe design
    with quality filtering.
    
    Args:
        gene_name (str): Name of the target gene
        transcriptome (Transcriptome): Transcriptome object with gene annotation data
        species_identifier (str): Species identifier for database organization
        base_dir (str, optional): Base directory for input files and output. Defaults to "".
        
    Returns:
        int: Number of available high-quality probe pairs
        
    Raises:
        ValueError: If gene is not found in transcriptome
        Exception: If any step of the pipeline fails
        
    Notes:
        - Uses longest CDS transcript for probe design
        - Performs BLAST analysis against both mature and pre-mRNA databases
        - Masks off-target binding regions
        - Designs probes using B1 amplifier by default
        - Prints progress information during execution
        
    Examples:
        >>> n_probes = check_probe_availability("gene1", transcriptome, "dmel", "/data")
        Number of available probes for gene1: 45
        >>> print(f"Can design {n_probes} probe pairs")
    """
    try:
        # Retrieve gene from transcriptome
        gene = transcriptome.get_gene(gene_name)
        if gene is None:
            raise ValueError(f"Gene '{gene_name}' not found in transcriptome")
        
        # Get the target sequence 
        if not hasattr(gene, 'target_sequence') or gene.target_sequence is None:
            raise ValueError(f"Gene '{gene_name}' does not have a target_sequence attribute")

        # Perform BLAST analysis and process results to identify off-target regions
        print(f"Running BLAST analysis for {gene_name}...")
        blast_gene(
            gene_name,
            transcriptome,
            species_identifier=species_identifier,
            base_dir=base_dir,
            permitted_off_targets=permitted_off_targets
        )

        # Design probes on unique regions using B1 amplifier
        print(f"Designing HCR-FISH probes for {gene_name}...")
        probes, regions, positions = design_hcr_probes(gene.unique_sequence, amplifier="B1")

        # Report results
        n_probes = len(probes)
        print(f"Number of available probes for {gene_name}: {n_probes}")
        
        if n_probes == 0:
            print(f"Warning: No suitable probe regions found for {gene_name}")
            print("This may be due to:")
            print("  - Extensive off-target binding sites")
            print("  - Poor sequence quality (gaps, extreme GC content)")
            print("  - Short target sequence")
        
        return n_probes

    except ValueError as e:
        print(f"Error: {e}")
        raise
    except Exception as e:
        print(f"Unexpected error during probe availability check for {gene_name}: {e}")
        raise Exception(f"Probe availability check failed for {gene_name}: {e}")
    


def assign_target(
        gene_name: str,
        transcriptome: Transcriptome,
        sequence_type: str = 'mRNA',
        transcript_id: Optional[str] = None,
        use_longest_cds: bool = True,
        use_longest_bounds: bool = False
) -> None:
    """
    Assign a target sequence to a gene from its transcriptome.
    
    This function selects a specific transcript for a given gene and extracts 
    either its mRNA or DNA sequence as the target for probe design. The transcript 
    can be selected by ID, or automatically as the longest CDS or bounds.
    
    Args:
        gene_name (str): Name or identifier of the target gene
        transcriptome (Transcriptome): Transcriptome object containing gene annotations
        sequence_type (str, optional): Type of sequence to extract. Either 'mRNA' 
            (mature transcript without introns) or 'DNA' (genomic DNA with introns). 
            Defaults to 'mRNA'.
        transcript_id (Optional[str], optional): Specific transcript ID to use. If provided,
            this transcript will be selected regardless of use_longest_cds or 
            use_longest_bounds settings. Defaults to None.
        use_longest_cds (bool, optional): If True, select the transcript with the 
            longest coding sequence (CDS). Cannot be used with use_longest_bounds or 
            transcript_id. Defaults to True.
        use_longest_bounds (bool, optional): If True, select the transcript spanning 
            the longest genomic region. Cannot be used with use_longest_cds or 
            transcript_id. Defaults to False.
    
    Returns:
        None: Stores the selected sequence in gene.target_sequence attribute.
    
    Raises:
        ValueError: If sequence_type is not 'mRNA' or 'DNA'
        ValueError: If both use_longest_cds and use_longest_bounds are True
        ValueError: If transcript_id is specified along with use_longest_cds or use_longest_bounds
        ValueError: If gene is not found in transcriptome
        ValueError: If gene has no transcripts available
        ValueError: If no transcript selection method is specified
        ValueError: If specified transcript_id is not found for the gene
        
    Examples:
        >>> # Use longest coding sequence (default)
        >>> assign_target("Or9a", transcriptome)
        Assigned mRNA sequence for gene 'Or9a' from transcript 'NM_001001234.1' with length 2341 bp.
        
        >>> # Use longest genomic span
        >>> assign_target("dsx", transcriptome, use_longest_cds=False, use_longest_bounds=True)
        Assigned mRNA sequence for gene 'dsx' from transcript 'NM_001001235.1' with length 2156 bp.
        
        >>> # Extract genomic DNA sequence
        >>> assign_target("Or9a", transcriptome, sequence_type='DNA', use_longest_cds=True)
        Assigned DNA sequence for gene 'Or9a' from transcript 'NM_001001234.1' with length 2500 bp.
        
        >>> # Use specific transcript by ID
        >>> assign_target("Or9a", transcriptome, transcript_id="NM_001001234.2")
        Assigned mRNA sequence for gene 'Or9a' from transcript 'NM_001001234.2' with length 1950 bp.
    
    Notes:
        - The selected transcript sequence is stored in gene.target_sequence for use by 
          downstream functions like design_hcr_probes() or blast_gene().
        - Exactly one of use_longest_cds, use_longest_bounds, or transcript_id must be specified.
        - mRNA sequences are typically preferred for probe design as they represent the 
          mature transcript that will be targeted.
        - DNA sequences include introns and may be useful for designing probes that span 
          exon-intron junctions.
        - The function prints informational output showing the selected transcript and 
          sequence length.
    """
    # Validate sequence type
    if sequence_type not in ['mRNA', 'DNA']:
        raise ValueError("Invalid sequence_type. Must be 'mRNA' or 'DNA'.")

    # Validate mutually exclusive parameters
    if use_longest_cds and use_longest_bounds:
        raise ValueError("Cannot use both longest CDS and longest bounds simultaneously.")
    
    if transcript_id and (use_longest_cds or use_longest_bounds):
        raise ValueError("Cannot specify transcript_id when using longest CDS or bounds.")
    
    # Validate gene exists
    gene = transcriptome.get_gene(gene_name)
    if gene is None:
        raise ValueError(f"Gene '{gene_name}' not found in transcriptome.")
    
    # Validate gene has transcripts
    if len(gene.transcripts) == 0:
        raise ValueError(f"Gene '{gene_name}' has no transcripts available.")
    
    # Validate at least one selection method is provided
    if not (use_longest_cds or use_longest_bounds or transcript_id):
        raise ValueError("Must specify one of use_longest_cds, use_longest_bounds, or transcript_id.")
    
    # Select transcript based on specified method
    if use_longest_cds:
        transcript = gene.get_transcript_longest_cds()
    elif use_longest_bounds:
        transcript = gene.get_transcript_longest_bounds()
    elif transcript_id:
        transcript = gene.get_transcript(transcript_id)
        if transcript is None:
            raise ValueError(f"Transcript ID '{transcript_id}' not found for gene '{gene_name}'.")
    
    # Extract appropriate sequence type
    sequence = transcript.mrna_sequence if sequence_type == 'mRNA' else transcript.dna_sequence

    # Assign to gene object
    gene.target_sequence = sequence
    print(f"Assigned {sequence_type} sequence for gene '{gene_name}' from transcript '{transcript.name}' with length {len(sequence)} bp.")

    
