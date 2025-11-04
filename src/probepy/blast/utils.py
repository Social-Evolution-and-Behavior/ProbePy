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
import logging
from typing import Union, Any, Optional
from probepy.blast.install import ensure_blast_tools

logger = logging.getLogger(__name__)

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
        logger.error(f"Input FASTA file not found: {input_fasta}")
        return False
    
    if dbtype not in ['nucl', 'prot']:
        logger.error(f"Invalid dbtype: {dbtype}. Must be 'nucl' or 'prot'")
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
        logger.info(f"Successfully created BLAST database: {output_name}")
        if result.stdout.strip():
            logger.debug(f"  {result.stdout.strip()}")
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error creating BLAST database:")
        logger.error(f"Command: {' '.join(cmd)}")
        logger.error(f"Return code: {e.returncode}")
        if e.stderr:
            logger.error(f"Error: {e.stderr.strip()}")
        return False
        
    except subprocess.TimeoutExpired:
        logger.error(f"makeblastdb timed out after 5 minutes")
        return False


def run_blastn(
    query: str,
    database: str,
    output_file: Optional[str] = None,
    evalue: float = 0.001,
    max_target_seqs: int = 10000,
    outfmt: str = "10 qseqid sseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    num_threads: int = 4,
    task: str = "blastn",
    ungapped: bool = False,
    word_size: Optional[int] = None,
    strand: str = "both",
    reward: Optional[int] = None,
    penalty: Optional[int] = None,
    dust: str = "yes",
    soft_masking: bool = True,
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
            to report. Defaults to 10000.
        outfmt (str, optional): Output format specification. Defaults to 
            "10 qseqid sseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore"
            (CSV format with specific columns).
        num_threads (int, optional): Number of threads to use for search.
            Defaults to 4.
        task (str, optional): BLAST task preset. Defaults to "blastn".
        ungapped (bool, optional): Perform ungapped alignment only. Defaults to False.
        word_size (Optional[int], optional): Word size for initial matches. Defaults to None (BLAST default).
        strand (str, optional): Query strand(s) to search. Options: "both", "plus", "minus". 
            Defaults to "both".
        reward (Optional[int], optional): Reward for nucleotide match. Defaults to None (BLAST default).
        penalty (Optional[int], optional): Penalty for nucleotide mismatch. Defaults to None (BLAST default).
        dust (str, optional): Filter query sequence with DUST. Options: "yes", "no". Defaults to "yes".
        soft_masking (bool, optional): Apply filtering locations as soft masks. Defaults to True.
        **kwargs (Any): Additional BLASTN parameters as key-value pairs,
            passed directly to the blastn command.
    
    Returns:
        Union[str, bool]: If output_file is None, returns search results as
            string. If output_file is specified, returns True on success.
            Returns False on any error.
                         
    Examples:
        >>> # Save results to file
        >>> success = run_blastn("probes.fasta", "transcriptome_db", "results.txt")
        
        >>> # Get results as string with custom parameters
        >>> results = run_blastn(
        ...     "probes.fasta", "transcriptome_db", 
        ...     ungapped=True, word_size=15, strand="plus"
        ... )
        >>> if results:
        ...     print(f"Found {len(results.splitlines())} hits")
    """
    if not ensure_blast_tools():
        return False
    
    # Validate inputs
    if not os.path.exists(query):
        logger.error(f"Query file not found: {query}")
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
        logger.error(f"BLAST database not found: {database}")
        logger.error(f"  Expected files like: {database}.nin, {database}.nhr, {database}.nsq")
        return False
    
    # Build command with common parameters
    cmd = [
        'blastn',
        '-query',
        query,
        '-db',
        database,
        '-task',
        task,
        '-evalue',
        str(evalue),
        '-max_target_seqs',
        str(max_target_seqs),
        '-outfmt',
        str(outfmt),
        '-num_threads',
        str(num_threads),
        '-strand',
        strand,
        '-dust',
        dust,
        '-soft_masking',
        str(soft_masking).lower()
    ]
    
    # Add optional parameters
    if ungapped:
        cmd.append('-ungapped')
    
    if word_size is not None:
        cmd.extend(['-word_size', str(word_size)])
    
    if reward is not None:
        cmd.extend(['-reward', str(reward)])
    
    if penalty is not None:
        cmd.extend(['-penalty', str(penalty)])
    
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
            logger.info(f"[OK] BLAST results saved to: {output_file}")
            return True
        else:
            return result.stdout
            
    except subprocess.CalledProcessError as e:
        logger.error(f"BLASTN search failed:")
        logger.error(f"Command: {' '.join(cmd)}")
        logger.error(f"Return code: {e.returncode}")
        if e.stderr:
            logger.error(f"Error: {e.stderr.strip()}")
        return False
        
    except subprocess.TimeoutExpired:
        logger.error(f"BLASTN search timed out after 10 minutes")
        logger.error("Consider reducing the query size or using fewer threads")
        return False