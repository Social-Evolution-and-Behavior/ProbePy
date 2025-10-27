#!/usr/bin/env python3
"""
Functions to automatically install NCBI command line tools (BLAST+) based on the operating system.

This module provides automated installation functionality for BLAST+ tools based on the
operating system. It ensures the tools are available for makeblastdb and blastn operations
in the HCR-FISH probe design pipeline.
"""

import os
import sys
import platform
import subprocess
import urllib.request
import tarfile
import zipfile
import shutil
from typing import Dict, Union, Any
from hcrfish.blast.config import BLAST_VERSION


def check_blast_tools() -> Dict[str, Dict[str, Union[bool, str, None]]]:
    """
    Check if BLAST+ tools are available and return their status.
    
    Iterates through required BLAST+ tools (makeblastdb and blastn) and checks
    if they are available in the system PATH. Returns version information for
    available tools and None for those not found.
    
    Returns:
        Dict[str, Dict[str, Union[bool, str, None]]]: Status dictionary with tool
            names as keys. Each tool has 'available' (bool) and 'version' (str or None)
            fields indicating installation status and version information.
            
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
                [
                    tool,
                    '-version'
                ],
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
            print(f"[OK] {tool}: {version_line}")
            
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            status[tool] = {
                'available': False,
                'version': None
            }
            print(f"[NOT FOUND] {tool}: Not found")
    
    return status


def ensure_blast_tools() -> bool:
    """
    Ensure BLAST+ tools are available. If not, attempt automatic installation.
    
    This function performs a quick check to see if BLAST+ tools are available
    in the system PATH. If the critical makeblastdb tool is not found, it notifies
    the user to run the installation function. This is a lightweight check that
    doesn't require full installation attempts.
    
    Returns:
        bool: True if tools are available and accessible, False otherwise.
        
    Examples:
        >>> if ensure_blast_tools():
        ...     print("BLAST+ tools ready to use")
        ... else:
        ...     print("Please install BLAST+ tools manually")
    """
    try:
        subprocess.run(
            [
                'makeblastdb',
                '-version'
            ],
            capture_output=True,
            check=True,
            timeout=10
        )
        return True
        
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        print("BLAST+ tools not found. Please run hcrfish.install_blast_tools() to install them.")
        return False


def get_system_info() -> str:
    """
    Get system information to determine the appropriate BLAST+ package.
    
    Detects the current operating system and architecture to identify which
    BLAST+ binary package should be downloaded and installed.
    
    Returns:
        str: System identifier string (e.g., 'macosx-universal', 'linux-x64', 'win64').
        
    Raises:
        RuntimeError: If the operating system is not supported.
        
    Examples:
        >>> system = get_system_info()
        >>> print(f"Detected system: {system}")
    """
    system = platform.system()
    machine = platform.machine()
    
    if system == "Darwin":  # macOS
        if machine == "arm64":
            return "macosx-universal"
        else:
            return "macosx-universal"
    elif system == "Linux":
        if machine == "x86_64":
            return "linux-x64"
        elif machine == "aarch64":
            return "linux-aarch64"
        else:
            return "linux-x64"  # fallback
    elif system == "Windows":
        return "win64"
    else:
        raise RuntimeError(f"Unsupported operating system: {system}")


def check_if_installed() -> bool:
    """
    Check if BLAST+ tools are already installed and accessible.
    
    Verifies that makeblastdb can be executed from the system PATH and reports
    the version information if installation is detected.
    
    Returns:
        bool: True if BLAST+ tools are installed and accessible, False otherwise.
        
    Examples:
        >>> if check_if_installed():
        ...     print("BLAST+ is ready to use")
    """
    try:
        result = subprocess.run(
            [
                "makeblastdb",
                "-version"
            ],
            capture_output=True,
            text=True,
            check=True
        )
        print("[OK] BLAST+ tools are already installed:")
        print(result.stdout.strip())
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False



def install_with_homebrew() -> bool:
    """
    Try to install BLAST+ using Homebrew on macOS.
    
    Attempts to detect Homebrew installation and use it to install the BLAST+
    package. If Homebrew is not available or installation fails, returns False
    to allow fallback to alternative installation methods.
    
    Returns:
        bool: True if installation was successful, False otherwise.
        
    Examples:
        >>> if install_with_homebrew():
        ...     print("BLAST+ installed with Homebrew")
    """
    try:
        # Check if homebrew is available
        subprocess.run(
            [
                "brew",
                "--version"
            ],
            capture_output=True,
            check=True
        )
        print("Installing BLAST+ using Homebrew...")
        
        result = subprocess.run(
            [
                "brew",
                "install",
                "blast"
            ],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print("Successfully installed BLAST+ using Homebrew")
            return True
        else:
            print(f"Homebrew installation failed: {result.stderr}")
            return False
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Homebrew not available, trying manual installation...")
        return False


def install_with_apt() -> bool:
    """
    Try to install BLAST+ using apt on Ubuntu/Debian systems.
    
    Attempts to use the system apt package manager to install ncbi-blast+.
    Requires sudo privileges. If apt is not available or installation fails,
    returns False to allow fallback to alternative installation methods.
    
    Returns:
        bool: True if installation was successful, False otherwise.
        
    Examples:
        >>> if install_with_apt():
        ...     print("BLAST+ installed with apt")
    """
    try:
        print("Installing BLAST+ using apt...")
        
        # Update package list
        subprocess.run(
            [
                "sudo",
                "apt",
                "update"
            ],
            check=True
        )
        
        # Install ncbi-blast+
        result = subprocess.run(
            [
                "sudo",
                "apt",
                "install",
                "-y",
                "ncbi-blast+"
            ],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print("Successfully installed BLAST+ using apt")
            return True
        else:
            print(f"Apt installation failed: {result.stderr}")
            return False
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Apt not available, trying manual installation...")
        return False


def manual_blast_install() -> bool:
    """
    Manually download and install BLAST+ from NCBI.
    
    Downloads BLAST+ binaries directly from the NCBI FTP server and extracts
    them to the user's local bin directory. This function is used as a fallback
    when system package managers are not available. For Windows, it provides
    instructions for manual installation.
    
    Returns:
        bool: True if installation was successful, False otherwise.
        
    Examples:
        >>> if manual_blast_install():
        ...     print("BLAST+ installed successfully")
    """
    print("Downloading BLAST+ directly from NCBI...")
    
    try:
        system_type = get_system_info()
        
        # Determine download URL and file extension
        if system_type.startswith("win"):
            filename = f"ncbi-blast-{BLAST_VERSION}+-{system_type}.exe"
            url = f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{BLAST_VERSION}/{filename}"
        else:
            filename = f"ncbi-blast-{BLAST_VERSION}+-{system_type}.tar.gz"
            url = f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{BLAST_VERSION}/{filename}"
        
        print(f"Downloading {filename}...")
        
        # Create temporary directory
        temp_dir = "temp_blast_install"
        os.makedirs(
            temp_dir,
            exist_ok=True
        )
        
        # Download file
        file_path = os.path.join(
            temp_dir,
            filename
        )
        urllib.request.urlretrieve(
            url,
            file_path
        )
        
        # Extract and install
        if system_type.startswith("win"):
            print("Windows installer downloaded. Please run the .exe file manually.")
            print(f"File location: {os.path.abspath(file_path)}")
            return False
        else:
            print("Extracting BLAST+ tools...")
            with tarfile.open(
                file_path,
                'r:gz'
            ) as tar:
                tar.extractall(temp_dir)
            
            # Find the extracted directory
            extracted_dir = None
            for item in os.listdir(temp_dir):
                item_path = os.path.join(
                    temp_dir,
                    item
                )
                if os.path.isdir(item_path) and item.startswith("ncbi-blast"):
                    extracted_dir = item_path
                    break
            
            if not extracted_dir:
                raise RuntimeError("Could not find extracted BLAST+ directory")
            
            # Install to user's local bin directory
            home_dir = os.path.expanduser("~")
            bin_dir = os.path.join(
                home_dir,
                ".local",
                "bin"
            )
            os.makedirs(
                bin_dir,
                parents=True,
                exist_ok=True
            )
            
            blast_bin_dir = os.path.join(
                extracted_dir,
                "bin"
            )
            for tool in os.listdir(blast_bin_dir):
                tool_path = os.path.join(
                    blast_bin_dir,
                    tool
                )
                if os.path.isfile(tool_path):
                    dest_path = os.path.join(
                        bin_dir,
                        tool
                    )
                    shutil.copy2(
                        tool_path,
                        dest_path
                    )
                    os.chmod(
                        dest_path,
                        0o755
                    )

            print(f"BLAST+ tools installed to {bin_dir}")
            print(f"Make sure {bin_dir} is in your PATH environment variable")

            # Add to PATH suggestion
            shell_rc = os.path.join(
                home_dir,
                ".bashrc"
            )
            if not os.path.exists(shell_rc):
                shell_rc = os.path.join(
                    home_dir,
                    ".zshrc"
                )
            
            shell_name = os.path.basename(shell_rc)
            print(f"\nTo add to PATH, add this line to your {shell_name}:")
            print(f'export PATH="$HOME/.local/bin:$PATH"')
        
        # Cleanup
        shutil.rmtree(temp_dir)
        return True
        
    except Exception as e:
        print(f"Manual installation failed: {e}")
        return False


def install_blast_tools() -> None:
    """
    Main installation function for BLAST+ tools.
    
    Orchestrates the complete installation process for BLAST+ tools. Checks if
    the tools are already installed, then attempts system-specific installation
    methods (Homebrew on macOS, apt on Linux) before falling back to manual
    installation from NCBI. Performs a final verification after installation
    completes.
    
    Returns:
        None: Prints status messages to stdout and stderr.
        
    Examples:
        >>> install_blast_tools()
        # Installs BLAST+ and provides status updates
    """
    
    # Check if already installed
    if check_if_installed():
        print("BLAST+ tools are already installed.")
        return

    print("BLAST+ tools not found. Installing...")

    system = platform.system()
    
    # Try different installation methods based on the system
    success = False
    
    # System-specific package managers
    if system == "Darwin":  # macOS
        if install_with_homebrew():
            success = True
    elif system == "Linux":
        if install_with_apt():
            success = True
    
    # Fallback to manual installation
    if not success:
        print("Trying manual installation...")
        success = manual_blast_install()
    
    # Final check
    if success:
        print("\nChecking installation...")
        if check_if_installed():
            print("BLAST+ tools successfully installed!")
        else:
            print("Installation completed but tools not found in PATH.")
            print("You may need to restart your terminal or update your PATH.")
    else:
        print("Installation failed. Please install BLAST+ manually from:")
        print("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download")

