"""
Utility functions for ensuring NCBI BLAST+ tools are available.

This module provides functions to check, install, and run NCBI BLAST+ tools
required for probe specificity verification in the HCR-FISH pipeline.
"""

import subprocess
import sys
import os
from pathlib import Path
from typing import Dict, Union, Any, Optional


def check_blast_tools() -> Dict[str, Dict[str, Union[bool, str, None]]]:
    """
    Check if BLAST+ tools are available and return their status.
    
    Checks for the presence and version of essential BLAST+ tools required
    for the HCR probe design pipeline.
    
    Returns:
        Dict[str, Dict[str, Union[bool, str, None]]]: Status dictionary with tool names as keys.
            Each tool has 'available' (bool) and 'version' (str or None) fields.
            
    Examples:
        >>> status = check_blast_tools()
        >>> if status['blastn']['available']:
        ...     print(f"BLASTN version: {status['blastn']['version']}")
        ... else:
        ...     print("BLASTN not found")
    """
    required_tools = ['makeblastdb', 'blastn']
    status: Dict[str, Dict[str, Union[bool, str, None]]] = {}
    
    for tool in required_tools:
        try:
            result = subprocess.run(
                [tool, '-version'], 
                capture_output=True, 
                text=True, 
                check=True,
                timeout=30
            )
            version_line = result.stdout.strip().split('\n')[0]
            status[tool] = {
                'available': True,
                'version': version_line
            }
            print(f"✓ {tool}: {version_line}")
            
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            status[tool] = {
                'available': False,
                'version': None
            }
            print(f"❌ {tool}: Not found")
    
    return status


def ensure_blast_tools() -> bool:
    """
    Ensure BLAST+ tools are available. If not, attempt automatic installation.
    
    This function first checks if BLAST+ tools are available in the system PATH.
    If not found, it attempts to run the included installation script to download
    and install the tools automatically.
    
    Returns:
        bool: True if tools are available after check/installation, False otherwise
        
    Examples:
        >>> if ensure_blast_tools():
        ...     print("BLAST+ tools ready to use")
        ... else:
        ...     print("Please install BLAST+ tools manually")
    """
    # Quick check for most critical tool
    try:
        subprocess.run(
            ['makeblastdb', '-version'], 
            capture_output=True, 
            check=True,
            timeout=10
        )
        return True
        
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        print("BLAST+ tools not found. Attempting automatic installation...")
        
        # Locate installation script
        script_dir = Path(__file__).parent
        install_script = script_dir / "install_ncbi_tools.py"
        
        if not install_script.exists():
            print(f"❌ Installation script not found at {install_script}")
            print("Please install BLAST+ tools manually:")
            print("  macOS: brew install blast")
            print("  Ubuntu/Debian: sudo apt-get install ncbi-blast+")
            print("  Windows: Download from https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/")
            return False
        
        try:
            # Run installation script
            print(f"Running installation script: {install_script}")
            result = subprocess.run(
                [sys.executable, str(install_script)],
                timeout=300  # 5 minute timeout for installation
            )
            
            if result.returncode != 0:
                print("❌ Installation script failed")
                return False
                
        except subprocess.TimeoutExpired:
            print("❌ Installation timed out")
            return False
        except Exception as e:
            print(f"❌ Installation error: {e}")
            return False
        
        # Verify installation worked
        try:
            subprocess.run(
                ['makeblastdb', '-version'], 
                capture_output=True, 
                check=True,
                timeout=10
            )
            print("✓ BLAST+ tools successfully installed and verified")
            return True
            
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            print("❌ Installation completed but tools still not accessible")
            print("You may need to restart your kernel/terminal or add BLAST+ to your PATH")
            return False


def run_makeblastdb(input_fasta: Union[str, Path], 
                   output_name: Union[str, Path],
                   dbtype: str = 'nucl', 
                   parse_seqids: bool = True,
                   title: Optional[str] = None) -> bool:
    """
    Create a BLAST database from a FASTA file.
    
    Args:
        input_fasta (Union[str, Path]): Path to input FASTA file
        output_name (Union[str, Path]): Output database name/path
        dbtype (str, optional): Database type ('nucl' or 'prot'). Defaults to 'nucl'.
        parse_seqids (bool, optional): Whether to parse sequence IDs. Defaults to True.
        title (Optional[str], optional): Database title for documentation. Defaults to None.
    
    Returns:
        bool: True if database creation successful, False otherwise
        
    Examples:
        >>> success = run_makeblastdb("transcripts.fasta", "transcripts_db")
        >>> if success:
        ...     print("Database ready for searching")
    """
    if not ensure_blast_tools():
        return False
    
    # Validate inputs
    input_path = Path(input_fasta)
    if not input_path.exists():
        print(f"❌ Input FASTA file not found: {input_path}")
        return False
    
    if dbtype not in ['nucl', 'prot']:
        print(f"❌ Invalid dbtype: {dbtype}. Must be 'nucl' or 'prot'")
        return False
    
    # Build command
    cmd = ['makeblastdb', '-in', str(input_path), '-dbtype', dbtype, '-out', str(output_name)]
    
    if parse_seqids:
        cmd.append('-parse_seqids')
        
    if title:
        cmd.extend(['-title', title])
    
    try:
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            check=True,
            timeout=300  # 5 minute timeout
        )
        print(f"✓ Successfully created BLAST database: {output_name}")
        if result.stdout.strip():
            print(f"  {result.stdout.strip()}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error creating BLAST database:")
        print(f"  Command: {' '.join(cmd)}")
        print(f"  Return code: {e.returncode}")
        if e.stderr:
            print(f"  Error: {e.stderr.strip()}")
        return False
        
    except subprocess.TimeoutExpired:
        print(f"❌ makeblastdb timed out after 5 minutes")
        return False


def run_blastn(query: Union[str, Path], 
               database: Union[str, Path],
               output_file: Optional[Union[str, Path]] = None,
               evalue: float = 0.001,
               max_target_seqs: int = 500,
               outfmt: str = "6",
               num_threads: int = 1,
               **kwargs: Any) -> Union[str, bool]:
    """
    Run BLASTN search with comprehensive error handling.
    
    Args:
        query (Union[str, Path]): Path to query sequences (FASTA format)
        database (Union[str, Path]): Path to BLAST database (without extension)
        output_file (Optional[Union[str, Path]], optional): Output file path. If None, returns results as string.
        evalue (float, optional): E-value threshold. Defaults to 0.001.
        max_target_seqs (int, optional): Maximum number of target sequences. Defaults to 500.
        outfmt (str, optional): Output format. Defaults to "6" (tabular).
        num_threads (int, optional): Number of threads to use. Defaults to 1.
        **kwargs (Any): Additional BLASTN parameters as key-value pairs.
    
    Returns:
        Union[str, bool]: If output_file is None, returns search results as string.
                         If output_file is specified, returns True on success.
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
    query_path = Path(query)
    if not query_path.exists():
        print(f"❌ Query file not found: {query_path}")
        return False
    
    # Check if database files exist
    db_path = Path(database)
    db_extensions = ['.nin', '.nhr', '.nsq']  # Required nucleotide database files
    
    if not any((db_path.parent / (db_path.name + ext)).exists() for ext in db_extensions):
        print(f"❌ BLAST database not found: {database}")
        print(f"  Expected files like: {database}.nin, {database}.nhr, {database}.nsq")
        return False
    
    # Build command with common parameters
    cmd = [
        'blastn',
        '-query', str(query_path),
        '-db', str(database),
        '-evalue', str(evalue),
        '-max_target_seqs', str(max_target_seqs),
        '-outfmt', str(outfmt),
        '-num_threads', str(num_threads)
    ]
    
    # Add output file if specified
    if output_file:
        cmd.extend(['-out', str(output_file)])
    
    # Add additional parameters from kwargs
    for key, value in kwargs.items():
        if key.startswith('-'):
            cmd.extend([key, str(value)])
        else:
            cmd.extend([f'-{key}', str(value)])
    
    try:
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            check=True,
            timeout=600  # 10 minute timeout
        )
        
        if output_file:
            print(f"✓ BLAST results saved to: {output_file}")
            return True
        else:
            return result.stdout
            
    except subprocess.CalledProcessError as e:
        print(f"❌ BLASTN search failed:")
        print(f"  Command: {' '.join(cmd)}")
        print(f"  Return code: {e.returncode}")
        if e.stderr:
            print(f"  Error: {e.stderr.strip()}")
        return False
        
    except subprocess.TimeoutExpired:
        print(f"❌ BLASTN search timed out after 10 minutes")
        print("  Consider reducing the query size or using fewer threads")
        return False