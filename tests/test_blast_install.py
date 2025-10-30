"""
Tests for probepy.blast.install module.

This module tests BLAST+ installation and tool checking functionality.
"""

import pytest
import subprocess
import tempfile
import os
import platform
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open

from probepy.blast.install import (
    check_blast_tools,
    ensure_blast_tools,
    get_system_info,
    check_if_installed,
    install_with_homebrew,
    install_with_apt,
    manual_blast_install,
    install_blast_tools
)


class TestCheckBlastTools:
    """Test BLAST tool availability checking."""
    
    @patch('subprocess.run')
    def test_tools_available(self, mock_run):
        """Test when BLAST tools are available."""
        # Mock successful version check
        mock_result = MagicMock()
        mock_result.stdout = "makeblastdb: 2.15.0+\nPackage: blast 2.15.0"
        mock_result.returncode = 0
        mock_run.return_value = mock_result
        
        status = check_blast_tools()
        
        assert isinstance(status, dict)
        assert 'makeblastdb' in status
        assert 'blastn' in status
        assert status['makeblastdb']['available'] is True
        assert status['blastn']['available'] is True
        assert 'version' in status['makeblastdb']
        assert 'makeblastdb' in status['makeblastdb']['version']
    
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
        assert status['blastn']['available'] is False
        assert status['blastn']['version'] is None
    
    @patch('subprocess.run')
    def test_timeout_handling(self, mock_run):
        """Test timeout handling in tool checking."""
        mock_run.side_effect = subprocess.TimeoutExpired(['makeblastdb', '-version'], 30)
        
        status = check_blast_tools()
        
        assert status['makeblastdb']['available'] is False
        assert status['makeblastdb']['version'] is None
    
    @patch('subprocess.run')
    def test_called_process_error_handling(self, mock_run):
        """Test handling of CalledProcessError."""
        mock_run.side_effect = subprocess.CalledProcessError(1, ['makeblastdb', '-version'])
        
        status = check_blast_tools()
        
        assert status['makeblastdb']['available'] is False
        assert status['makeblastdb']['version'] is None
    
    @patch('subprocess.run')
    def test_partial_tool_availability(self, mock_run):
        """Test when only one tool is available."""
        def mock_run_side_effect(*args, **kwargs):
            if 'makeblastdb' in str(args[0]):
                result = MagicMock()
                result.stdout = "makeblastdb: 2.15.0+"
                result.returncode = 0
                return result
            else:
                raise FileNotFoundError()
        
        mock_run.side_effect = mock_run_side_effect
        
        status = check_blast_tools()
        
        assert status['makeblastdb']['available'] is True
        assert status['blastn']['available'] is False


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
    def test_tools_not_available(self, mock_run):
        """Test when tools are not available."""
        # Mock command not found
        mock_run.side_effect = FileNotFoundError()
        
        result = ensure_blast_tools()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_timeout_during_check(self, mock_run):
        """Test timeout during tool check."""
        mock_run.side_effect = subprocess.TimeoutExpired(['makeblastdb', '-version'], 10)
        
        result = ensure_blast_tools()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_called_process_error_during_check(self, mock_run):
        """Test CalledProcessError during tool check."""
        mock_run.side_effect = subprocess.CalledProcessError(1, ['makeblastdb', '-version'])
        
        result = ensure_blast_tools()
        
        assert result is False


class TestGetSystemInfo:
    """Test system information detection."""
    
    @patch('platform.system')
    @patch('platform.machine')
    def test_macos_arm64(self, mock_machine, mock_system):
        """Test macOS ARM64 (Apple Silicon) detection."""
        mock_system.return_value = "Darwin"
        mock_machine.return_value = "arm64"
        
        result = get_system_info()
        
        assert result == "macosx-universal"
    
    @patch('platform.system')
    @patch('platform.machine')
    def test_macos_x86_64(self, mock_machine, mock_system):
        """Test macOS x86_64 (Intel) detection."""
        mock_system.return_value = "Darwin"
        mock_machine.return_value = "x86_64"
        
        result = get_system_info()
        
        assert result == "macosx-universal"
    
    @patch('platform.system')
    @patch('platform.machine')
    def test_linux_x86_64(self, mock_machine, mock_system):
        """Test Linux x86_64 detection."""
        mock_system.return_value = "Linux"
        mock_machine.return_value = "x86_64"
        
        result = get_system_info()
        
        assert result == "linux-x64"
    
    @patch('platform.system')
    @patch('platform.machine')
    def test_linux_aarch64(self, mock_machine, mock_system):
        """Test Linux ARM64 detection."""
        mock_system.return_value = "Linux"
        mock_machine.return_value = "aarch64"
        
        result = get_system_info()
        
        assert result == "linux-aarch64"
    
    @patch('platform.system')
    @patch('platform.machine')
    def test_linux_fallback(self, mock_machine, mock_system):
        """Test Linux fallback for unknown architecture."""
        mock_system.return_value = "Linux"
        mock_machine.return_value = "ppc64le"
        
        result = get_system_info()
        
        assert result == "linux-x64"  # Fallback
    
    @patch('platform.system')
    @patch('platform.machine')
    def test_windows(self, mock_machine, mock_system):
        """Test Windows detection."""
        mock_system.return_value = "Windows"
        mock_machine.return_value = "AMD64"
        
        result = get_system_info()
        
        assert result == "win64"
    
    @patch('platform.system')
    def test_unsupported_os(self, mock_system):
        """Test unsupported operating system raises error."""
        mock_system.return_value = "FreeBSD"
        
        with pytest.raises(RuntimeError, match="Unsupported operating system"):
            get_system_info()


class TestCheckIfInstalled:
    """Test checking if BLAST tools are installed."""
    
    @patch('subprocess.run')
    def test_already_installed(self, mock_run):
        """Test when BLAST is already installed."""
        mock_result = MagicMock()
        mock_result.stdout = "makeblastdb: 2.15.0+"
        mock_result.returncode = 0
        mock_run.return_value = mock_result
        
        result = check_if_installed()
        
        assert result is True
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_not_installed(self, mock_run):
        """Test when BLAST is not installed."""
        mock_run.side_effect = FileNotFoundError()
        
        result = check_if_installed()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_called_process_error(self, mock_run):
        """Test CalledProcessError handling."""
        mock_run.side_effect = subprocess.CalledProcessError(1, ['makeblastdb', '-version'])
        
        result = check_if_installed()
        
        assert result is False


class TestInstallWithHomebrew:
    """Test Homebrew installation method."""
    
    @patch('subprocess.run')
    def test_homebrew_not_available(self, mock_run):
        """Test when Homebrew is not available."""
        mock_run.side_effect = FileNotFoundError()
        
        result = install_with_homebrew()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_successful_installation(self, mock_run):
        """Test successful installation with Homebrew."""
        # First call checks brew version, second installs blast
        def mock_run_side_effect(*args, **kwargs):
            if 'brew' in str(args[0]) and '--version' in str(args[0]):
                return MagicMock(returncode=0, stdout="Homebrew 4.0.0")
            else:
                return MagicMock(returncode=0, stdout="Installing blast...")
        
        mock_run.side_effect = mock_run_side_effect
        
        result = install_with_homebrew()
        
        assert result is True
        assert mock_run.call_count == 2
    
    @patch('subprocess.run')
    def test_installation_failure(self, mock_run):
        """Test failed installation with Homebrew."""
        # First call checks brew version, second fails
        def mock_run_side_effect(*args, **kwargs):
            if 'brew' in str(args[0]) and '--version' in str(args[0]):
                return MagicMock(returncode=0, stdout="Homebrew 4.0.0")
            else:
                return MagicMock(returncode=1, stderr="Installation failed", stdout="")
        
        mock_run.side_effect = mock_run_side_effect
        
        result = install_with_homebrew()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_brew_check_error(self, mock_run):
        """Test error during brew version check."""
        mock_run.side_effect = subprocess.CalledProcessError(1, ['brew', '--version'])
        
        result = install_with_homebrew()
        
        assert result is False


class TestInstallWithApt:
    """Test apt installation method."""
    
    @patch('subprocess.run')
    def test_apt_not_available(self, mock_run):
        """Test when apt is not available."""
        mock_run.side_effect = FileNotFoundError()
        
        result = install_with_apt()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_successful_installation(self, mock_run):
        """Test successful installation with apt."""
        # All calls succeed
        mock_run.return_value = MagicMock(returncode=0, stdout="Installing...")
        
        result = install_with_apt()
        
        assert result is True
        assert mock_run.call_count == 2  # apt update + apt install
    
    @patch('subprocess.run')
    def test_installation_failure(self, mock_run):
        """Test failed installation with apt."""
        # Update succeeds, install fails
        def mock_run_side_effect(*args, **kwargs):
            if 'update' in str(args[0]):
                return MagicMock(returncode=0)
            else:
                return MagicMock(returncode=1, stderr="Install failed", stdout="")
        
        mock_run.side_effect = mock_run_side_effect
        
        result = install_with_apt()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_update_failure(self, mock_run):
        """Test failure during apt update."""
        mock_run.side_effect = subprocess.CalledProcessError(1, ['sudo', 'apt', 'update'])
        
        result = install_with_apt()
        
        assert result is False


class TestManualBlastInstall:
    """Test manual BLAST installation."""
    
    @patch('probepy.blast.install.get_system_info')
    @patch('urllib.request.urlretrieve')
    @patch('os.makedirs')
    @patch('tarfile.open')
    @patch('shutil.rmtree')
    @patch('shutil.copy2')
    @patch('os.chmod')
    @patch('os.listdir')
    @patch('os.path.isdir')
    @patch('os.path.isfile')
    @patch('os.path.exists')
    @patch('os.path.expanduser')
    def test_successful_linux_installation(self, mock_expanduser, mock_exists, mock_isfile, 
                                          mock_isdir, mock_listdir, mock_chmod, mock_copy,
                                          mock_rmtree, mock_tarfile, mock_makedirs,
                                          mock_urlretrieve, mock_system_info):
        """Test successful manual installation on Linux."""
        mock_system_info.return_value = "linux-x64"
        mock_expanduser.return_value = "/home/user"
        mock_exists.return_value = True
        mock_isdir.return_value = True
        mock_isfile.return_value = True
        mock_listdir.side_effect = [
            ["ncbi-blast-2.15.0+"],  # temp_dir contents
            ["makeblastdb", "blastn"]  # bin directory contents
        ]
        
        # Mock tarfile extraction
        mock_tar = MagicMock()
        mock_tarfile.return_value.__enter__.return_value = mock_tar
        
        result = manual_blast_install()
        
        assert result is True
        mock_urlretrieve.assert_called_once()
        mock_tarfile.assert_called_once()
        assert mock_copy.call_count == 2  # Copy 2 tools
    
    @patch('probepy.blast.install.get_system_info')
    def test_windows_installation_message(self, mock_system_info):
        """Test Windows installation returns False with message."""
        mock_system_info.return_value = "win64"
        
        with patch('urllib.request.urlretrieve'):
            with patch('os.makedirs'):
                with patch('os.path.abspath', return_value="/path/to/installer.exe"):
                    result = manual_blast_install()
        
        assert result is False  # Windows requires manual installation
    
    @patch('probepy.blast.install.get_system_info')
    @patch('urllib.request.urlretrieve')
    def test_download_failure(self, mock_urlretrieve, mock_system_info):
        """Test failed download."""
        mock_system_info.return_value = "linux-x64"
        mock_urlretrieve.side_effect = Exception("Download failed")
        
        with patch('os.makedirs'):
            result = manual_blast_install()
        
        assert result is False
    
    @patch('probepy.blast.install.get_system_info')
    @patch('urllib.request.urlretrieve')
    @patch('os.makedirs')
    @patch('tarfile.open')
    def test_extraction_failure(self, mock_tarfile, mock_makedirs, 
                               mock_urlretrieve, mock_system_info):
        """Test failed tarfile extraction."""
        mock_system_info.return_value = "linux-x64"
        mock_tarfile.side_effect = Exception("Extraction failed")
        
        result = manual_blast_install()
        
        assert result is False


class TestInstallBlastTools:
    """Test main installation orchestration function."""
    
    @patch('probepy.blast.install.check_if_installed')
    def test_already_installed(self, mock_check):
        """Test when BLAST is already installed."""
        mock_check.return_value = True
        
        install_blast_tools()
        
        # Should return early without attempting installation
        mock_check.assert_called_once()
    
    @patch('probepy.blast.install.check_if_installed')
    @patch('probepy.blast.install.install_with_homebrew')
    @patch('platform.system')
    def test_macos_homebrew_success(self, mock_platform, mock_homebrew, mock_check):
        """Test successful installation on macOS with Homebrew."""
        mock_check.side_effect = [False, True]  # Not installed, then installed
        mock_platform.return_value = "Darwin"
        mock_homebrew.return_value = True
        
        install_blast_tools()
        
        mock_homebrew.assert_called_once()
    
    @patch('probepy.blast.install.check_if_installed')
    @patch('probepy.blast.install.install_with_apt')
    @patch('platform.system')
    def test_linux_apt_success(self, mock_platform, mock_apt, mock_check):
        """Test successful installation on Linux with apt."""
        mock_check.side_effect = [False, True]  # Not installed, then installed
        mock_platform.return_value = "Linux"
        mock_apt.return_value = True
        
        install_blast_tools()
        
        mock_apt.assert_called_once()
    
    @patch('probepy.blast.install.check_if_installed')
    @patch('probepy.blast.install.install_with_homebrew')
    @patch('probepy.blast.install.manual_blast_install')
    @patch('platform.system')
    def test_fallback_to_manual(self, mock_platform, mock_manual, mock_homebrew, mock_check):
        """Test fallback to manual installation when package manager fails."""
        mock_check.side_effect = [False, True]  # Not installed, then installed
        mock_platform.return_value = "Darwin"
        mock_homebrew.return_value = False  # Homebrew fails
        mock_manual.return_value = True  # Manual succeeds
        
        install_blast_tools()
        
        mock_homebrew.assert_called_once()
        mock_manual.assert_called_once()
    
    @patch('probepy.blast.install.check_if_installed')
    @patch('probepy.blast.install.install_with_homebrew')
    @patch('probepy.blast.install.manual_blast_install')
    @patch('platform.system')
    def test_all_installation_methods_fail(self, mock_platform, mock_manual, 
                                          mock_homebrew, mock_check):
        """Test when all installation methods fail."""
        mock_check.return_value = False  # Never installed
        mock_platform.return_value = "Darwin"
        mock_homebrew.return_value = False
        mock_manual.return_value = False
        
        install_blast_tools()
        
        # Should try both methods
        mock_homebrew.assert_called_once()
        mock_manual.assert_called_once()
    
    @patch('probepy.blast.install.check_if_installed')
    @patch('probepy.blast.install.install_with_apt')
    @patch('probepy.blast.install.manual_blast_install')
    @patch('platform.system')
    def test_installation_success_but_not_in_path(self, mock_platform, mock_manual,
                                                  mock_apt, mock_check):
        """Test when installation succeeds but tools not found in PATH."""
        mock_check.return_value = False  # Not found in PATH
        mock_platform.return_value = "Linux"
        mock_apt.return_value = False
        mock_manual.return_value = True  # Installation reports success
        
        install_blast_tools()
        
        # Should attempt installation and final check
        assert mock_check.call_count >= 2


class TestIntegration:
    """Integration tests for installation workflow."""
    
    @patch('subprocess.run')
    def test_check_and_ensure_workflow(self, mock_run):
        """Test complete check and ensure workflow."""
        # First check fails
        mock_run.side_effect = FileNotFoundError()
        
        # Check tools
        status = check_blast_tools()
        assert all(not tool['available'] for tool in status.values())
        
        # Ensure tools (should fail)
        result = ensure_blast_tools()
        assert result is False
    
    @patch('subprocess.run')
    def test_full_installation_simulation(self, mock_run):
        """Test full installation simulation."""
        # Simulate: not installed -> install -> installed
        call_count = [0]
        
        def mock_run_side_effect(*args, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                # First check: not installed
                raise FileNotFoundError()
            else:
                # After "installation": installed
                return MagicMock(returncode=0, stdout="makeblastdb: 2.15.0+")
        
        mock_run.side_effect = mock_run_side_effect
        
        # Initial check
        assert not ensure_blast_tools()
        
        # After "installation"
        assert ensure_blast_tools()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
