"""
Tests for probepy.hcr.prep module.

This module tests the preparation functions including rsync downloads,
mRNA export to FASTA, and BLAST database creation.
"""

import os
import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import subprocess

from probepy.hcr.prep import download_with_rsync, export_mrna_to_fasta, create_blast_databases
from probepy.transcriptomics.classes import Transcriptome, Gene, Transcript


class TestDownloadWithRsync:
    """Test download_with_rsync function."""
    
    def test_invalid_file_type_raises_error(self):
        """Test that invalid file_type raises ValueError."""
        with pytest.raises(ValueError, match="file_type must be either 'genome' or 'transcriptome'"):
            download_with_rsync("rsync://example.com/file.fa", "test_species", "invalid_type")
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    @patch('pathlib.Path.mkdir')
    def test_successful_download_genome(self, mock_mkdir, mock_exists, mock_subprocess):
        """Test successful genome download."""
        # Setup mocks - file doesn't exist initially, then exists after download
        mock_exists.side_effect = [False, True]  # First check: file doesn't exist, second: file exists after download
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            result = download_with_rsync(
                "rsync://example.com/test.fa",
                "test_species",
                "genome",
                base_dir=temp_dir
            )
            
            # Verify rsync was called - check first call (rsync should be first)
            mock_subprocess.assert_called()
            calls = mock_subprocess.call_args_list
            rsync_call = calls[0][0][0]  # First call's first argument
            assert "rsync" in rsync_call
            assert "--partial" in rsync_call
            assert "--progress" in rsync_call
            assert "rsync://example.com/test.fa" in rsync_call
            
            # Verify result path
            expected_path = os.path.join(temp_dir, "input", "test_species", "genome", "test.fa")
            assert result == expected_path
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    @patch('pathlib.Path.mkdir')
    def test_successful_download_with_decompression(self, mock_mkdir, mock_exists, mock_subprocess):
        """Test successful download with .gz decompression."""
        # Setup mocks - file doesn't exist initially, .gz exists after download, .fa exists after gunzip
        # Checks: final_path.exists (False), output_path.exists (True after rsync), decompressed.exists (False), then True after gunzip
        mock_exists.side_effect = [False, True, False, True]
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            result = download_with_rsync(
                "rsync://example.com/test.fa.gz",
                "test_species",
                "genome",
                base_dir=temp_dir
            )
            
            # Verify both rsync and gunzip were called
            assert mock_subprocess.call_count >= 2  # rsync, gunzip (chflags might not be called on non-macOS)
            
            # Verify result path (without .gz extension)
            expected_path = os.path.join(temp_dir, "input", "test_species", "genome", "test.fa")
            assert result == expected_path
    
    @patch('subprocess.run')
    def test_rsync_failure_raises_error(self, mock_subprocess):
        """Test that rsync failure raises CalledProcessError."""
        mock_subprocess.side_effect = subprocess.CalledProcessError(
            1, ["rsync"], "rsync failed"
        )
        
        with tempfile.TemporaryDirectory() as temp_dir:
            with pytest.raises(subprocess.CalledProcessError):
                download_with_rsync(
                    "rsync://example.com/test.fa",
                    "test_species",
                    "genome",
                    base_dir=temp_dir
                )
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    @patch('pathlib.Path.mkdir')
    def test_missing_file_after_download_raises_error(self, mock_mkdir, mock_exists, mock_subprocess):
        """Test that missing file after download raises FileNotFoundError."""
        # rsync succeeds but file is not found
        mock_exists.return_value = False
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            with pytest.raises(FileNotFoundError, match="Downloaded file not found"):
                download_with_rsync(
                    "rsync://example.com/test.fa",
                    "test_species",
                    "genome",
                    base_dir=temp_dir
                )
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    @patch('pathlib.Path.mkdir')
    @patch('pathlib.Path.rename')
    def test_fna_to_fa_conversion(self, mock_rename, mock_mkdir, mock_exists, mock_subprocess):
        """Test .fna to .fa file extension conversion."""
        # File doesn't exist initially, .fna exists after download, .fa doesn't exist yet for rename
        # Checks: final_path.exists (False), output_path.exists (True after rsync), new_path.exists (False for rename)
        mock_exists.side_effect = [False, True, False]
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            result = download_with_rsync(
                "rsync://example.com/test.fna",
                "test_species",
                "genome",
                base_dir=temp_dir
            )
            
            # Verify rename was called
            mock_rename.assert_called_once()
            
            # Verify result has .fa extension
            assert result.endswith(".fa")
            assert "test.fa" in result


class TestExportMrnaToFasta:
    """Test export_mrna_to_fasta function."""
    
    def test_empty_transcriptome_raises_error(self):
        """Test that empty transcriptome raises ValueError."""
        empty_transcriptome = Transcriptome()
        
        with pytest.raises(ValueError, match="Transcriptome object is empty"):
            export_mrna_to_fasta(empty_transcriptome, "test_species")
    
    def test_successful_export(self):
        """Test successful mRNA export to FASTA files."""
        # Create test transcriptome
        transcriptome = Transcriptome()
        
        # Create test gene with transcript
        gene = Gene("FBgn0000001", "test_gene")
        gene.chromosome = "chr2L"
        gene.strand = "+"
        
        transcript = Transcript("NM_001")
        transcript.chromosome = "chr2L"
        transcript.strand = "+"
        transcript.position = (1000, 2000)
        transcript.mrna_sequence = "ATCGATCGATCG"
        transcript.dna_sequence = "ATCGATCGATCGATCGATCG"  # includes introns
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            no_introns_path, yes_introns_path = export_mrna_to_fasta(
                transcriptome, "test_species", base_dir=temp_dir
            )
            
            # Verify files were created
            assert Path(no_introns_path).exists()
            assert Path(yes_introns_path).exists()
            
            # Check file contents
            with open(no_introns_path, 'r') as f:
                content = f.read()
                assert "test_gene" in content
                assert "ATCGATCGATCG" in content
                assert "chr2L:1000-2000" in content
                assert "strand=+" in content
            
            with open(yes_introns_path, 'r') as f:
                content = f.read()
                assert "test_gene" in content
                assert "ATCGATCGATCGATCGATCG" in content
    
    def test_no_sequences_raises_error(self):
        """Test that transcriptome with no sequences raises ValueError."""
        transcriptome = Transcriptome()
        
        # Create gene with transcript but no sequences
        gene = Gene("FBgn0000001", "test_gene")
        transcript = Transcript("NM_001")
        transcript.mrna_sequence = ""
        transcript.dna_sequence = ""
        
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            with pytest.raises(ValueError, match="No transcripts with sequence data found"):
                export_mrna_to_fasta(transcriptome, "test_species", base_dir=temp_dir)
    
    def test_multiple_genes_and_transcripts(self):
        """Test export with multiple genes and transcripts."""
        transcriptome = Transcriptome()
        
        # Create multiple genes with transcripts
        for i in range(3):
            gene = Gene(f"FBgn000000{i}", f"test_gene_{i}")
            
            for j in range(2):
                transcript = Transcript(f"NM_00{i}_{j}")
                transcript.chromosome = f"chr{i}"
                transcript.strand = "+" if i % 2 == 0 else "-"
                transcript.position = (1000 + i*100, 2000 + i*100)
                transcript.mrna_sequence = "ATCG" * (10 + i)
                transcript.dna_sequence = "ATCG" * (15 + i)
                
                gene.add_transcript(transcript)
            
            transcriptome.add_gene(gene)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            no_introns_path, yes_introns_path = export_mrna_to_fasta(
                transcriptome, "test_species", base_dir=temp_dir
            )
            
            # Count sequences in files
            with open(no_introns_path, 'r') as f:
                no_introns_content = f.read()
                assert no_introns_content.count('>') == 6  # 3 genes × 2 transcripts
            
            with open(yes_introns_path, 'r') as f:
                yes_introns_content = f.read()
                assert yes_introns_content.count('>') == 6


class TestCreateBlastDatabases:
    """Test create_blast_databases function."""
    
    def test_missing_fasta_files_raises_error(self):
        """Test that missing FASTA files raise FileNotFoundError."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Don't create the FASTA files
            with pytest.raises(FileNotFoundError, match="FASTA file not found"):
                create_blast_databases(
                    species_identifier="test_species",
                    base_dir=temp_dir
                )
    
    @patch('subprocess.run')
    def test_makeblastdb_not_available_raises_error(self, mock_subprocess):
        """Test that missing makeblastdb raises FileNotFoundError."""
        mock_subprocess.side_effect = FileNotFoundError("makeblastdb not found")
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create the expected directory structure and FASTA files
            no_introns_dir = Path(temp_dir) / "input" / "test_species" / "transcriptome" / "mRNA_no_introns"
            yes_introns_dir = Path(temp_dir) / "input" / "test_species" / "transcriptome" / "mRNA_yes_introns"
            no_introns_dir.mkdir(parents=True, exist_ok=True)
            yes_introns_dir.mkdir(parents=True, exist_ok=True)
            
            fasta1 = no_introns_dir / "mRNA_no_introns.fasta"
            fasta2 = yes_introns_dir / "mRNA_yes_introns.fasta"
            fasta1.write_text(">seq1\nATCG")
            fasta2.write_text(">seq2\nGCTA")
            
            with pytest.raises(FileNotFoundError, match="makeblastdb not found"):
                create_blast_databases(species_identifier="test_species", base_dir=temp_dir)
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    def test_successful_database_creation(self, mock_exists, mock_subprocess):
        """Test successful BLAST database creation."""
        # Mock makeblastdb version check and database creation
        mock_subprocess.return_value = MagicMock(
            returncode=0,
            stdout="makeblastdb: 2.13.0+"
        )
        # Mock database files exist after creation
        mock_exists.return_value = True
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create the expected directory structure and FASTA files
            no_introns_dir = Path(temp_dir) / "input" / "test_species" / "transcriptome" / "mRNA_no_introns"
            yes_introns_dir = Path(temp_dir) / "input" / "test_species" / "transcriptome" / "mRNA_yes_introns"
            no_introns_dir.mkdir(parents=True, exist_ok=True)
            yes_introns_dir.mkdir(parents=True, exist_ok=True)
            
            fasta1 = no_introns_dir / "mRNA_no_introns.fasta"
            fasta2 = yes_introns_dir / "mRNA_yes_introns.fasta"
            fasta1.write_text(">seq1\nATCG")
            fasta2.write_text(">seq2\nGCTA")
            
            create_blast_databases(
                "test_species",
                base_dir=temp_dir
            )
            
            # Verify makeblastdb was called for both databases
            assert mock_subprocess.call_count == 3  # version + 2 databases
    
    @patch('subprocess.run')
    def test_makeblastdb_failure_raises_error(self, mock_subprocess):
        """Test that makeblastdb failure raises CalledProcessError."""
        # Version check succeeds, but database creation fails
        mock_subprocess.side_effect = [
            MagicMock(returncode=0, stdout="makeblastdb: 2.13.0+"),  # version check
            subprocess.CalledProcessError(1, ["makeblastdb"], "makeblastdb failed"),  # first db creation
        ]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create the expected directory structure and FASTA files
            no_introns_dir = Path(temp_dir) / "input" / "test_species" / "transcriptome" / "mRNA_no_introns"
            yes_introns_dir = Path(temp_dir) / "input" / "test_species" / "transcriptome" / "mRNA_yes_introns"
            no_introns_dir.mkdir(parents=True, exist_ok=True)
            yes_introns_dir.mkdir(parents=True, exist_ok=True)
            
            fasta1 = no_introns_dir / "mRNA_no_introns.fasta"
            fasta2 = yes_introns_dir / "mRNA_yes_introns.fasta"
            fasta1.write_text(">seq1\nATCG")
            fasta2.write_text(">seq2\nGCTA")
            
            with pytest.raises(subprocess.CalledProcessError):
                create_blast_databases(species_identifier="test_species", base_dir=temp_dir)


class TestIntegrationWorkflow:
    """Integration tests for the complete preparation workflow."""
    
    @patch('subprocess.run')
    def test_complete_workflow(self, mock_subprocess):
        """Test the complete workflow from download to database creation."""
        # Setup mocks
        mock_subprocess.return_value = MagicMock(returncode=0, stdout="makeblastdb: 2.13.0+")
        
        # Create test transcriptome
        transcriptome = Transcriptome()
        gene = Gene("FBgn0000001", "test_gene")
        transcript = Transcript("NM_001")
        transcript.chromosome = "chr2L"
        transcript.strand = "+"
        transcript.position = (1000, 2000)
        transcript.mrna_sequence = "ATCGATCGATCG"
        transcript.dna_sequence = "ATCGATCGATCGATCGATCG"
        gene.add_transcript(transcript)
        transcriptome.add_gene(gene)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Step 1: Download (simulated with mocked rsync and file existence)
            with patch('pathlib.Path.exists', return_value=True):
                genome_path = download_with_rsync(
                    "rsync://example.com/test.fa.gz",
                    "test_species",
                    "genome",
                    base_dir=temp_dir
                )
            
            # Step 2: Export mRNA to FASTA (use real directories - don't mock mkdir)
            fasta_paths = export_mrna_to_fasta(transcriptome, "test_species", base_dir=temp_dir)
            
            # Step 3: Create BLAST databases (mocked)
            with patch('pathlib.Path.exists', return_value=True):
                create_blast_databases("test_species", base_dir=temp_dir)
            
            # Verify all steps completed
            assert genome_path is not None
            assert len(fasta_paths) == 2
            
            # Verify subprocess calls were made
            assert mock_subprocess.call_count >= 3  # rsync, gunzip, chflags, version, 2 makeblastdb


@pytest.fixture
def sample_transcriptome():
    """Create a sample transcriptome for testing."""
    transcriptome = Transcriptome()
    
    # Add a gene with transcript
    gene = Gene("FBgn0000001", "Or9a")
    gene.chromosome = "chr2L"
    gene.strand = "+"
    
    transcript = Transcript("NM_001001234.1")
    transcript.chromosome = "chr2L"
    transcript.strand = "+"
    transcript.position = (1000, 2000)
    transcript.mrna_sequence = ("ATGAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
                               "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                               "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
                               "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATA")
    transcript.dna_sequence = transcript.mrna_sequence + "GTAAGTAAGTAAGTAAGTAA"  # add introns
    transcript.cds_length = 194
    
    gene.add_transcript(transcript)
    transcriptome.add_gene(gene)
    
    return transcriptome


class TestWithSampleData:
    """Tests using sample transcriptome data."""
    
    def test_export_with_sample_transcriptome(self, sample_transcriptome):
        """Test export function with sample transcriptome."""
        with tempfile.TemporaryDirectory() as temp_dir:
            no_introns_path, yes_introns_path = export_mrna_to_fasta(
                sample_transcriptome, "test_species", base_dir=temp_dir
            )
            
            # Verify files exist and contain expected content
            assert Path(no_introns_path).exists()
            assert Path(yes_introns_path).exists()
            
            with open(no_introns_path, 'r') as f:
                content = f.read()
                assert "Or9a" in content
                assert "NM_001001234.1" in content
                assert "chr2L:1000-2000" in content