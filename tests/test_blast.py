"""
Tests for hcrfish.blast module.

This module tests BLAST-related functionality including tool checking,
database creation, and sequence searching.
"""

import pytest
import subprocess
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, MagicMock

from hcrfish.blast.blast_utils import (
    check_blast_tools,
    ensure_blast_tools, 
    run_makeblastdb,
    run_blastn
)


class TestCheckBlastTools:
    """Test BLAST tool availability checking."""
    
    @patch('subprocess.run')
    def test_tools_available(self, mock_run):
        """Test when BLAST tools are available."""
        # Mock successful version check
        mock_result = MagicMock()
        mock_result.stdout = "makeblastdb: 2.12.0+\n Build: Mar 8 2022"
        mock_result.returncode = 0
        mock_run.return_value = mock_result
        
        status = check_blast_tools()
        
        assert isinstance(status, dict)
        assert 'makeblastdb' in status
        assert 'blastn' in status
        assert status['makeblastdb']['available'] is True
        assert 'version' in status['makeblastdb']
    
    @patch('subprocess.run')
    def test_tools_not_available(self, mock_run):
        """Test when BLAST tools are not available."""
        # Mock command not found
        mock_run.side_effect = FileNotFoundError()
        
        status = check_blast_tools()
        
        assert isinstance(status, dict)
        assert 'makeblastdb' in status
        assert 'blastn' in status
        assert status['makeblastdb']['available'] is False
        assert status['makeblastdb']['version'] is None
    
    @patch('subprocess.run')
    def test_timeout_handling(self, mock_run):
        """Test timeout handling in tool checking."""
        mock_run.side_effect = subprocess.TimeoutExpired(['makeblastdb', '-version'], 30)
        
        status = check_blast_tools()
        
        assert status['makeblastdb']['available'] is False
        assert status['makeblastdb']['version'] is None


class TestEnsureBlastTools:
    """Test BLAST tool installation and verification."""
    
    @patch('subprocess.run')
    def test_tools_already_available(self, mock_run):
        """Test when tools are already available."""
        # Mock successful version check
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_run.return_value = mock_result
        
        result = ensure_blast_tools()
        
        assert result is True
        # Should only call version check, not installation
        assert mock_run.call_count == 1
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    def test_installation_script_missing(self, mock_exists, mock_run):
        """Test when installation script is missing."""
        # First call fails (tools not available), exists check fails
        mock_run.side_effect = FileNotFoundError()
        mock_exists.return_value = False
        
        result = ensure_blast_tools()
        
        assert result is False
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    def test_installation_success(self, mock_exists, mock_run):
        """Test successful installation."""
        # Set up mock sequence: fail, script exists, install succeeds, verify succeeds
        mock_exists.return_value = True
        
        def mock_run_side_effect(*args, **kwargs):
            if 'makeblastdb' in str(args[0]):
                if mock_run.call_count <= 1:
                    raise FileNotFoundError()  # First check fails
                else:
                    return MagicMock(returncode=0)  # After install succeeds
            else:
                return MagicMock(returncode=0)  # Installation script succeeds
        
        mock_run.side_effect = mock_run_side_effect
        
        result = ensure_blast_tools()
        
        assert result is True
    
    @patch('subprocess.run')
    @patch('pathlib.Path.exists')
    def test_installation_timeout(self, mock_exists, mock_run):
        """Test installation timeout."""
        mock_exists.return_value = True
        
        def mock_run_side_effect(*args, **kwargs):
            if 'makeblastdb' in str(args[0]):
                raise FileNotFoundError()
            else:
                raise subprocess.TimeoutExpired(args[0], 300)
        
        mock_run.side_effect = mock_run_side_effect
        
        result = ensure_blast_tools()
        
        assert result is False


class TestRunMakeblastdb:
    """Test makeblastdb wrapper function."""
    
    def test_input_validation(self):
        """Test input file validation."""
        # Non-existent input file
        result = run_makeblastdb("/nonexistent/file.fasta", "test_db")
        assert result is False
    
    def test_dbtype_validation(self, temp_fasta_file):
        """Test database type validation."""
        result = run_makeblastdb(temp_fasta_file, "test_db", dbtype="invalid")
        assert result is False
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    def test_blast_tools_not_available(self, mock_ensure, temp_fasta_file):
        """Test when BLAST tools are not available."""
        mock_ensure.return_value = False
        
        result = run_makeblastdb(temp_fasta_file, "test_db")
        
        assert result is False
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('subprocess.run')
    def test_successful_database_creation(self, mock_run, mock_ensure, temp_fasta_file):
        """Test successful database creation."""
        mock_ensure.return_value = True
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Building a new DB"
        mock_run.return_value = mock_result
        
        result = run_makeblastdb(temp_fasta_file, "test_db")
        
        assert result is True
        mock_run.assert_called_once()
        
        # Check command structure
        call_args = mock_run.call_args[0][0]
        assert 'makeblastdb' in call_args
        assert '-in' in call_args
        assert '-dbtype' in call_args
        assert '-out' in call_args
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('subprocess.run')
    def test_database_creation_failure(self, mock_run, mock_ensure, temp_fasta_file):
        """Test database creation failure."""
        mock_ensure.return_value = True
        mock_run.side_effect = subprocess.CalledProcessError(1, ['makeblastdb'])
        
        result = run_makeblastdb(temp_fasta_file, "test_db")
        
        assert result is False
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('subprocess.run')
    def test_timeout_handling(self, mock_run, mock_ensure, temp_fasta_file):
        """Test timeout handling in database creation."""
        mock_ensure.return_value = True
        mock_run.side_effect = subprocess.TimeoutExpired(['makeblastdb'], 300)
        
        result = run_makeblastdb(temp_fasta_file, "test_db")
        
        assert result is False


class TestRunBlastn:
    """Test BLASTN wrapper function."""
    
    def test_query_file_validation(self):
        """Test query file validation."""
        result = run_blastn("/nonexistent/query.fasta", "database")
        assert result is False
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    def test_blast_tools_not_available(self, mock_ensure, temp_fasta_file):
        """Test when BLAST tools are not available."""
        mock_ensure.return_value = False
        
        result = run_blastn(temp_fasta_file, "database")
        
        assert result is False
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('pathlib.Path.exists')
    def test_database_not_found(self, mock_exists, mock_ensure, temp_fasta_file):
        """Test when database files don't exist."""
        mock_ensure.return_value = True
        mock_exists.return_value = False
        
        result = run_blastn(temp_fasta_file, "nonexistent_db")
        
        assert result is False
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('pathlib.Path.exists')
    @patch('subprocess.run')
    def test_successful_search_to_file(self, mock_run, mock_exists, mock_ensure, temp_fasta_file, temp_directory):
        """Test successful BLAST search with output file."""
        mock_ensure.return_value = True
        mock_exists.return_value = True
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "BLAST results"
        mock_run.return_value = mock_result
        
        output_file = os.path.join(temp_directory, "results.txt")
        result = run_blastn(temp_fasta_file, "test_db", output_file)
        
        assert result is True
        mock_run.assert_called_once()
        
        # Check command structure
        call_args = mock_run.call_args[0][0]
        assert 'blastn' in call_args
        assert '-query' in call_args
        assert '-db' in call_args
        assert '-out' in call_args
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('pathlib.Path.exists')
    @patch('subprocess.run')
    def test_successful_search_return_results(self, mock_run, mock_exists, mock_ensure, temp_fasta_file):
        """Test successful BLAST search returning results."""
        mock_ensure.return_value = True
        mock_exists.return_value = True
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "query1\tsubject1\t100.0\t50\t0\t0\t1\t50\t1\t50\t1e-20\t100"
        mock_run.return_value = mock_result
        
        result = run_blastn(temp_fasta_file, "test_db")
        
        assert isinstance(result, str)
        assert "query1" in result
        assert "subject1" in result
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('pathlib.Path.exists')
    @patch('subprocess.run')
    def test_custom_parameters(self, mock_run, mock_exists, mock_ensure, temp_fasta_file):
        """Test BLAST search with custom parameters."""
        mock_ensure.return_value = True
        mock_exists.return_value = True
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "results"
        mock_run.return_value = mock_result
        
        result = run_blastn(
            temp_fasta_file, 
            "test_db",
            evalue=1e-10,
            max_target_seqs=1000,
            word_size=11,
            reward=2,
            penalty=-3
        )
        
        assert isinstance(result, str)
        
        # Check that custom parameters were included
        call_args = mock_run.call_args[0][0]
        assert '-evalue' in call_args
        assert '1e-10' in call_args
        assert '-max_target_seqs' in call_args
        assert '1000' in call_args
        assert '-word_size' in call_args
        assert '11' in call_args
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('pathlib.Path.exists')
    @patch('subprocess.run')
    def test_search_failure(self, mock_run, mock_exists, mock_ensure, temp_fasta_file):
        """Test BLAST search failure."""
        mock_ensure.return_value = True
        mock_exists.return_value = True
        mock_run.side_effect = subprocess.CalledProcessError(1, ['blastn'])
        
        result = run_blastn(temp_fasta_file, "test_db")
        
        assert result is False
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('pathlib.Path.exists')
    @patch('subprocess.run')
    def test_search_timeout(self, mock_run, mock_exists, mock_ensure, temp_fasta_file):
        """Test BLAST search timeout."""
        mock_ensure.return_value = True
        mock_exists.return_value = True
        mock_run.side_effect = subprocess.TimeoutExpired(['blastn'], 600)
        
        result = run_blastn(temp_fasta_file, "test_db")
        
        assert result is False


class TestIntegration:
    """Integration tests for BLAST functionality."""
    
    @patch('hcrfish.blast.blast_utils.ensure_blast_tools')
    @patch('subprocess.run')
    def test_database_creation_and_search_pipeline(self, mock_run, mock_ensure, temp_fasta_file, temp_directory):
        """Test complete database creation and search pipeline."""
        mock_ensure.return_value = True
        
        # Mock successful database creation
        def mock_run_side_effect(*args, **kwargs):
            if 'makeblastdb' in str(args[0]):
                result = MagicMock()
                result.returncode = 0
                result.stdout = "Building a new DB"
                return result
            elif 'blastn' in str(args[0]):
                result = MagicMock()
                result.returncode = 0
                result.stdout = "query1\tsubject1\t100.0\t50"
                return result
            else:
                return MagicMock(returncode=0)
        
        mock_run.side_effect = mock_run_side_effect
        
        # Create database
        db_name = os.path.join(temp_directory, "test_db")
        db_success = run_makeblastdb(temp_fasta_file, db_name)
        assert db_success is True
        
        # Mock database files exist for search
        with patch('pathlib.Path.exists', return_value=True):
            # Search database
            search_result = run_blastn(temp_fasta_file, db_name)
            assert isinstance(search_result, str)
            assert "query1" in search_result
    
    def test_parameter_validation_chain(self, temp_fasta_file):
        """Test that parameter validation works correctly in sequence."""
        # Test invalid input file
        assert run_makeblastdb("/nonexistent.fasta", "db") is False
        
        # Test invalid database type
        assert run_makeblastdb(temp_fasta_file, "db", dbtype="invalid") is False
        
        # Test invalid query file for BLAST
        assert run_blastn("/nonexistent.fasta", "db") is False