"""
Preparation utilities for HCR-FISH probe design pipeline.

This module contains functions for downloading genomic data, exporting transcriptome
sequences to FASTA format, and creating BLAST databases for probe design workflows.
"""

import os
import subprocess
from typing import Optional
from pathlib import Path

from probepy.transcriptomics.classes import Transcriptome


def download_with_rsync(
    rsync_path: str,
    species_identifier: str,
    file_type: str = "genome",
    output_filename: Optional[str] = None,
    base_dir: str = "",
    overwrite: bool = False
):
    """
    Download genomic data files using rsync.
    
    This function downloads genome or transcriptome files from remote repositories
    using rsync, with automatic decompression and file organization by species.
    
    Args:
        rsync_path (str): Full rsync URL/path to the remote file
        species_identifier (str): Species identifier for directory organization (e.g., 'dmel', 'dyak')
        file_type (str, optional): Type of data ('genome' or 'transcriptome'). Defaults to 'genome'.
        output_filename (str, optional): Custom output filename. If None, infers from rsync_path.
        base_dir (str, optional): Base directory for downloads. Defaults to 'input'.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.
        
    Raises:
        ValueError: If file_type is not 'genome' or 'transcriptome'
        subprocess.CalledProcessError: If rsync or decompression fails
        FileNotFoundError: If the downloaded file is not found
        
    Examples:
        >>> # Download Drosophila melanogaster genome
        >>> genome_path = download_with_rsync(
        ...     "rsync://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz",
        ...     "dmel", 
        ...     "genome"
        ... )
        >>> 
        >>> # Download transcriptome annotation
        >>> gtf_path = download_with_rsync(
        ...     "rsync://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ncbiRefSeq.gtf.gz",
        ...     "dmel",
        ...     "transcriptome"
        ... )
        
    Notes:
        - Creates directory structure: {base_dir}/{species_identifier}/{file_type}/
        - Automatically decompresses .gz files
        - Converts .fna files to .fa extension for consistency
        - Uses rsync --partial --progress for robust downloads with progress reporting
        - On macOS, removes hidden file flags for better visibility
        - Skips download if final processed file already exists (unless overwrite=True)
    """
    if file_type not in ["genome", "transcriptome"]:
        raise ValueError("file_type must be either 'genome' or 'transcriptome'")
    
    # Create output directory
    output_dir = os.path.join(base_dir, "input", species_identifier, file_type)
    os.makedirs(output_dir, exist_ok=True)

    # Determine output filename
    if output_filename is None:
        output_filename = rsync_path.split('/')[-1]
    
    output_path = Path(output_dir) / output_filename

    # Determine the final processed file path (after decompression and renaming)
    final_path = Path(output_path)
    
    # Handle .gz decompression prediction
    if final_path.suffix == '.gz':
        final_path = final_path.with_suffix('')
    
    # Handle .fna to .fa conversion prediction
    if final_path.suffix == '.fna':
        final_path = final_path.with_suffix('.fa')
    
    # Check if final file already exists
    if final_path.exists() and not overwrite:
        print(f"File already exists: {final_path}")
        print("Skipping download (use overwrite=True to force re-download)")
        return str(final_path)
    
    print(f"Downloading {rsync_path}")
    print(f"Destination: {output_path}")
    if final_path != output_path:
        print(f"Final processed file will be: {final_path}")
    
    # Download file with rsync
    rsync_command = [
        "rsync",
        "--partial",
        "--progress",
        rsync_path,
        str(output_path)
    ]
    
    try:
        result = subprocess.run(
            rsync_command,
            check=True,
            capture_output=True,
            text=True
        )
        print("Download completed successfully")
        
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(
            e.returncode,
            e.cmd,
            f"rsync failed: {e.stderr}"
        )
    
    # Check if file was downloaded
    if not output_path.exists():
        raise FileNotFoundError(f"Downloaded file not found: {output_path}")
    
    # Handle decompression for .gz files
    if output_path.suffix == '.gz':
        decompressed_path = output_path.with_suffix('')
        
        # Check if decompressed file already exists
        if decompressed_path.exists() and not overwrite:
            print(f"Decompressed file already exists: {decompressed_path}")
            print("Removing downloaded .gz file")
            if output_path.exists():
                output_path.unlink()  # Remove the .gz file since we have the decompressed version
            output_path = decompressed_path
        else:
            print(f"Decompressing {output_path}")
            
            # If overwrite=True and decompressed file exists, remove it to avoid gunzip conflicts
            if overwrite and decompressed_path.exists():
                print(f"Removing existing decompressed file: {decompressed_path}")
                decompressed_path.unlink()
            
            # Check if the .gz file actually exists before trying to decompress
            if not output_path.exists():
                raise FileNotFoundError(f"Cannot decompress: .gz file not found at {output_path}")
            
            # Check if the file is actually gzipped by trying to read the magic number
            try:
                with open(output_path, 'rb') as f:
                    magic = f.read(2)
                    if magic != b'\x1f\x8b':
                        print(f"Warning: File {output_path} does not appear to be gzipped (magic number: {magic.hex()})")
                        print("Treating as already decompressed and renaming...")
                        # Rename the file to remove .gz extension
                        output_path.rename(decompressed_path)
                        output_path = decompressed_path
                        # Skip to the next processing step
                    else:
                        # File is properly gzipped, proceed with decompression
                        gunzip_command = ["gunzip", str(output_path)]
                        try:
                            result = subprocess.run(gunzip_command, check=True, capture_output=True, text=True)
                            # Update path to point to decompressed file
                            output_path = decompressed_path
                            print(f"Successfully decompressed to: {output_path}")
                            
                        except subprocess.CalledProcessError as e:
                            # Provide more detailed error information
                            error_msg = f"Decompression failed for {output_path}"
                            if e.stderr:
                                error_msg += f": {e.stderr.strip()}"
                            if e.stdout:
                                error_msg += f" (stdout: {e.stdout.strip()})"
                            
                            # Check if the file might already be decompressed
                            if decompressed_path.exists():
                                print(f"Warning: gunzip failed, but decompressed file exists: {decompressed_path}")
                                print("Using existing decompressed file")
                                if output_path.exists():
                                    output_path.unlink()  # Remove the problematic .gz file
                                output_path = decompressed_path
                            else:
                                raise subprocess.CalledProcessError(
                                    e.returncode,
                                    e.cmd,
                                    error_msg
                                )
            except Exception as e:
                print(f"Warning: Could not check file format: {e}")
                # Fall back to trying gunzip anyway
                gunzip_command = ["gunzip", str(output_path)]
                try:
                    result = subprocess.run(gunzip_command, check=True, capture_output=True, text=True)
                    output_path = decompressed_path
                    print(f"Successfully decompressed to: {output_path}")
                except subprocess.CalledProcessError as e:
                    error_msg = f"Decompression failed for {output_path}: {e.stderr if e.stderr else str(e)}"
                    if decompressed_path.exists():
                        print(f"Warning: gunzip failed, but decompressed file exists: {decompressed_path}")
                        print("Using existing decompressed file")
                        if output_path.exists():
                            output_path.unlink()
                        output_path = decompressed_path
                    else:
                        raise subprocess.CalledProcessError(e.returncode, e.cmd, error_msg)
    
    # Convert .fna to .fa for consistency
    if output_path.suffix == '.fna':
        new_path = output_path.with_suffix('.fa')
        
        # Check if .fa file already exists
        if new_path.exists() and not overwrite:
            print(f"Converted file already exists: {new_path}")
            print("Removing .fna file")
            output_path.unlink()  # Remove the .fna file since we have the .fa version
            output_path = new_path
        else:
            output_path.rename(new_path)
            output_path = new_path
            print(f"Renamed to {output_path}")
    
    # Remove hidden file flags on macOS
    try:
        subprocess.run(
            ["chflags", "nohidden", str(output_path)],
            check=False,  # Don't fail if chflags is not available
            capture_output=True
        )
    except FileNotFoundError:
        # chflags not available (not on macOS)
        pass
    
    print(f"File ready at: {output_path}")
    return str(output_path)


def export_mrna_to_fasta(
    transcriptome: Transcriptome,
    species_identifier: str,
    base_dir: str = "",
    overwrite: bool = False
):
    """
    Export mRNA sequences to FASTA files for BLAST database creation.
    
    This function creates two FASTA files: one with mature mRNA sequences (no introns)
    and one with pre-mRNA sequences (including introns). Both files are formatted
    for BLAST database creation and probe design workflows.
    
    Args:
        transcriptome (Transcriptome): Transcriptome object containing gene and transcript data
        species_identifier (str): Species identifier for directory organization (e.g., 'dmel', 'dyak')
        base_dir (str, optional): Base directory for output files. Defaults to 'input'.
        overwrite (bool, optional): Whether to overwrite existing FASTA files. Defaults to False.
        
    Raises:
        ValueError: If transcriptome is empty or contains no valid transcripts
        OSError: If output directories cannot be created or files cannot be written
        
    Examples:
        >>> from probepy.transcriptomics import load_transcriptome_object
        >>> transcriptome = load_transcriptome_object("dmel_transcriptome")
        >>> no_introns_path, yes_introns_path = export_mrna_to_fasta(transcriptome, "dmel")
        >>> print(f"Exported mature mRNA to: {no_introns_path}")
        >>> print(f"Exported pre-mRNA to: {yes_introns_path}")
        
    Notes:
        - Creates directory structure: {base_dir}/{species_identifier}/transcriptome/mRNA_no_introns/
        - Creates directory structure: {base_dir}/{species_identifier}/transcriptome/mRNA_yes_introns/
        - FASTA headers include transcript name, gene name, genomic location, and strand
        - Skips transcripts without sequence data
        - Files are compatible with makeblastdb for database creation
        - Skips export if files already exist (unless overwrite=True)
    """
    if not transcriptome.genes:
        raise ValueError("Transcriptome object is empty - no genes found")

    no_introns_dir = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_no_introns")
    yes_introns_dir = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_yes_introns")

    os.makedirs(no_introns_dir, exist_ok=True)
    os.makedirs(yes_introns_dir, exist_ok=True)
    
    # Define output file paths
    no_introns_path = os.path.join(no_introns_dir, "mRNA_no_introns.fasta")
    yes_introns_path = os.path.join(yes_introns_dir, "mRNA_yes_introns.fasta")
    
    # Check if files already exist and skip if overwrite=False
    if not overwrite:
        if Path(no_introns_path).exists() and Path(yes_introns_path).exists():
            print(f"FASTA files already exist for {species_identifier}")
            print("Skipping export (use overwrite=True to force re-export)")
            return (no_introns_path, yes_introns_path)

    print(f"Exporting mRNA sequences for {len(transcriptome.genes)} genes...")
    
    transcripts_no_introns = 0
    transcripts_yes_introns = 0
    
    try:
        # Export mRNA sequences without introns (mature mRNA)
        with open(no_introns_path, 'w') as no_introns_file:
            for gene_name, gene in transcriptome.genes.items():
                for transcript in gene.transcripts:
                    if transcript.mrna_sequence:
                        # Format: >transcript_id gene=gene_name location=chr:start-end strand=+/-
                        header = (
                            f">{transcript.name} gene={gene.name} "
                            f"location={transcript.chromosome}:{transcript.position[0]}-{transcript.position[1]} "
                            f"strand={transcript.strand}"
                        )
                        no_introns_file.write(f"{header}\n{transcript.mrna_sequence}\n")
                        transcripts_no_introns += 1
        
        # Export mRNA sequences with introns (pre-mRNA)
        with open(yes_introns_path, 'w') as yes_introns_file:
            for gene_name, gene in transcriptome.genes.items():
                for transcript in gene.transcripts:
                    if transcript.dna_sequence:
                        # Format: >transcript_id gene=gene_name location=chr:start-end strand=+/-
                        header = (
                            f">{transcript.name} gene={gene.name} "
                            f"location={transcript.chromosome}:{transcript.position[0]}-{transcript.position[1]} "
                            f"strand={transcript.strand}"
                        )
                        yes_introns_file.write(f"{header}\n{transcript.dna_sequence}\n")
                        transcripts_yes_introns += 1
                        
    except OSError as e:
        raise OSError(f"Failed to write FASTA files: {e}")
    
    print(f"Exported {transcripts_no_introns} transcripts to {no_introns_path}")
    print(f"Exported {transcripts_yes_introns} transcripts to {yes_introns_path}")
    
    if transcripts_no_introns == 0 and transcripts_yes_introns == 0:
        raise ValueError("No transcripts with sequence data found in transcriptome")
    
    return (no_introns_path, yes_introns_path)


def create_blast_databases(
    species_identifier: str,
    base_dir: str = ""
):
    """
    Create BLAST nucleotide databases from FASTA files.
    
    This function creates BLAST databases from mature mRNA and pre-mRNA FASTA files
    using the makeblastdb command. The databases are used for identifying potential
    off-target binding sites during probe design.
    
    Args:
        species_identifier (str): Species identifier for database organization (e.g., 'dmel', 'dyak')
        base_dir (str, optional): Base directory containing FASTA files. Defaults to ''.
        
    Returns:
        None: Databases are created in place; paths can be inferred from species_identifier
        
        
    Raises:
        FileNotFoundError: If FASTA files or makeblastdb command are not found
        subprocess.CalledProcessError: If makeblastdb commands fail
        ValueError: If FASTA files are empty or invalid
        
    Examples:
        >>> # After exporting FASTA files
        >>> db_paths = create_blast_databases("dmel", "/path/to/base")
        >>> no_introns_db, yes_introns_db = db_paths
        >>> print(f"BLAST databases created at: {db_paths}")
        
    Notes:
        - Requires BLAST+ tools to be installed and available in PATH
        - Creates databases with sequence ID parsing enabled (-parse_seqids)
        - Database files will have extensions (.ndb, .nhr, .nin, .not, .nsq, .ntf, .nto)
        - Output database names match the pattern used in blast_gene function
        - Verifies database creation by checking for required files
    """
    no_introns_fasta = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_no_introns", "mRNA_no_introns.fasta")
    yes_introns_fasta = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_yes_introns", "mRNA_yes_introns.fasta")

    # Verify input files exist
    if not Path(no_introns_fasta).exists():
        raise FileNotFoundError(f"FASTA file not found: {no_introns_fasta}")
    if not Path(yes_introns_fasta).exists():
        raise FileNotFoundError(f"FASTA file not found: {yes_introns_fasta}")
    
    # Check if makeblastdb is available
    try:
        result = subprocess.run(
            ["makeblastdb", "-version"],
            capture_output=True,
            text=True,
            check=True
        )
        print(f"Using makeblastdb version: {result.stdout.strip()}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise FileNotFoundError(
            "makeblastdb not found. Please install BLAST+ tools and ensure they are in your PATH."
        )
    
    # Define output database paths (without extensions)
    no_introns_db = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_no_introns", "mRNA_no_introns")
    yes_introns_db = os.path.join(base_dir, "input", species_identifier, "transcriptome", "mRNA_yes_introns", "mRNA_yes_introns")

    print("Creating BLAST databases...")
    
    # Create database for mature mRNA (no introns)
    print(f"Creating database: {no_introns_db}")
    no_introns_command = [
        "makeblastdb",
        "-in", no_introns_fasta,
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", no_introns_db
    ]
    
    try:
        result = subprocess.run(
            no_introns_command,
            capture_output=True,
            text=True,
            check=True
        )
        print("Mature mRNA database created successfully")
        
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(
            e.returncode,
            e.cmd,
            f"Failed to create no-introns database: {e.stderr}"
        )
    
    # Create database for pre-mRNA (with introns)
    print(f"Creating database: {yes_introns_db}")
    yes_introns_command = [
        "makeblastdb",
        "-in", yes_introns_fasta,
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", yes_introns_db
    ]
    
    try:
        result = subprocess.run(
            yes_introns_command,
            capture_output=True,
            text=True,
            check=True
        )
        print("Pre-mRNA database created successfully")
        
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(
            e.returncode,
            e.cmd,
            f"Failed to create introns database: {e.stderr}"
        )
    
    # Verify database files were created
    db_extensions = ['.ndb', '.nhr', '.nin', '.not', '.nsq', '.ntf', '.nto']
    
    for db_path in [no_introns_db, yes_introns_db]:
        missing_files = []
        for ext in db_extensions:
            if not Path(f"{db_path}{ext}").exists():
                missing_files.append(f"{db_path}{ext}")
        
        if missing_files:
            print(f"Warning: Some database files may be missing for {db_path}")
            print(f"Missing files: {missing_files}")
    
    print("BLAST database creation completed")
