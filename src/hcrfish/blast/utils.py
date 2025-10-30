"""
Utility functions for BLAST+ operations.

This module provides functions to run NCBI BLAST+ tools for probe specificity
verification in the HCR-FISH pipeline. It includes wrapper functions for
making BLAST databases and executing BLASTN searches with comprehensive
error handling and parameter management.
"""

import subprocess
import sys
import os
from typing import Union, Any, Optional
from probepy.blast.install import ensure_blast_tools

def run_makeblastdb(
    input_fasta: str,
    output_name: str,
    dbtype: str = 'nucl',
    parse_seqids: bool = True,
    title: Optional[str] = None
) -> bool:
    """
    Create a BLAST database from a FASTA file.
    
    Wraps the NCBI makeblastdb tool to create a searchable BLAST database
    from input FASTA sequences. Validates input parameters and provides
    comprehensive error reporting.
    
    Args:
        input_fasta (str): Path to input FASTA file.
        output_name (str): Output database name and path prefix.
        dbtype (str, optional): Database type, either 'nucl' for nucleotide
            or 'prot' for protein. Defaults to 'nucl'.
        parse_seqids (bool, optional): Whether to parse sequence identifiers
            in FASTA headers. Defaults to True.
        title (Optional[str], optional): Optional database title for
            documentation. Defaults to None.
    
    Returns:
        bool: True if database creation was successful, False otherwise.
        
    Examples:
        >>> success = run_makeblastdb("transcripts.fasta", "transcripts_db")
        >>> if success:
        ...     print("Database ready for searching")
    """
    if not ensure_blast_tools():
        return False
    
    # Validate inputs
    if not os.path.exists(input_fasta):
        print(f"Input FASTA file not found: {input_fasta}")
        return False
    
    if dbtype not in ['nucl', 'prot']:
        print(f"Invalid dbtype: {dbtype}. Must be 'nucl' or 'prot'")
        return False
    
    # Build command
    cmd = [
        'makeblastdb',
        '-in',
        input_fasta,
        '-dbtype',
        dbtype,
        '-out',
        output_name
    ]
    
    if parse_seqids:
        cmd.append('-parse_seqids')
        
    if title:
        cmd.extend([
            '-title',
            title
        ])
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=300  # 5 minute timeout
        )
        print(f"Successfully created BLAST database: {output_name}")
        if result.stdout.strip():
            print(f"  {result.stdout.strip()}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST database:")
        print(f"Command: {' '.join(cmd)}")
        print(f"Return code: {e.returncode}")
        if e.stderr:
            print(f"Error: {e.stderr.strip()}")
        return False
        
    except subprocess.TimeoutExpired:
        print(f"makeblastdb timed out after 5 minutes")
        return False


def run_blastn(
    query: str,
    database: str,
    output_file: Optional[str] = None,
    evalue: float = 0.001,
    max_target_seqs: int = 500,
    outfmt: str = "6",
    num_threads: int = 1,
    **kwargs: Any
) -> Union[str, bool]:
    """
    Run BLASTN search with comprehensive error handling.
    
    Wraps the NCBI blastn tool to search query sequences against a BLAST
    database. Can return results to a file or as a string. Includes extensive
    error handling and parameter validation.
    
    Args:
        query (str): Path to query sequences in FASTA format.
        database (str): Path to BLAST database without file extension.
        output_file (Optional[str], optional): Output file path. If None,
            returns results as string. Defaults to None.
        evalue (float, optional): E-value threshold for reporting hits.
            Defaults to 0.001.
        max_target_seqs (int, optional): Maximum number of target sequences
            to report. Defaults to 500.
        outfmt (str, optional): Output format code. Defaults to "6" (tabular).
        num_threads (int, optional): Number of threads to use for search.
            Defaults to 1.
        **kwargs (Any): Additional BLASTN parameters as key-value pairs,
            passed directly to the blastn command.
    
    Returns:
        Union[str, bool]: If output_file is None, returns search results as
            string. If output_file is specified, returns True on success.
            Returns False on any error.
                         
    Examples:
        >>> # Save results to file
        >>> success = run_blastn("probes.fasta", "transcriptome_db", "results.txt")
        
        >>> # Get results as string
        >>> results = run_blastn("probes.fasta", "transcriptome_db", outfmt="6")
        >>> if results:
        ...     print(f"Found {len(results.splitlines())} hits")
    """
    if not ensure_blast_tools():
        return False
    
    # Validate inputs
    if not os.path.exists(query):
        print(f"Query file not found: {query}")
        return False
    
    # Check if database files exist
    db_extensions = [
        '.nin',
        '.nhr',
        '.nsq'
    ]  # Required nucleotide database files
    
    db_found = False
    for ext in db_extensions:
        db_file = database + ext
        if os.path.exists(db_file):
            db_found = True
            break
    
    if not db_found:
        print(f"BLAST database not found: {database}")
        print(f"  Expected files like: {database}.nin, {database}.nhr, {database}.nsq")
        return False
    
    # Build command with common parameters
    cmd = [
        'blastn',
        '-query',
        query,
        '-db',
        database,
        '-evalue',
        str(evalue),
        '-max_target_seqs',
        str(max_target_seqs),
        '-outfmt',
        str(outfmt),
        '-num_threads',
        str(num_threads)
    ]
    
    # Add output file if specified
    if output_file:
        cmd.extend([
            '-out',
            output_file
        ])
    
    # Add additional parameters from kwargs
    for key, value in kwargs.items():
        if key.startswith('-'):
            cmd.extend([
                key,
                str(value)
            ])
        else:
            cmd.extend([
                f'-{key}',
                str(value)
            ])
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=600  # 10 minute timeout
        )
        
        if output_file:
            print(f"[OK] BLAST results saved to: {output_file}")
            return True
        else:
            return result.stdout
            
    except subprocess.CalledProcessError as e:
        print(f"BLASTN search failed:")
        print(f"Command: {' '.join(cmd)}")
        print(f"Return code: {e.returncode}")
        if e.stderr:
            print(f"Error: {e.stderr.strip()}")
        return False
        
    except subprocess.TimeoutExpired:
        print(f"BLASTN search timed out after 10 minutes")
        print("Consider reducing the query size or using fewer threads")
        return False