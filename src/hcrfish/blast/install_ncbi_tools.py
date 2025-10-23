#!/usr/bin/env python3
"""
Script to automatically install NCBI command line tools (BLAST+) based on the operating system.
This ensures the tools are available for makeblastdb and blastn operations.
"""

import os
import sys
import platform
import subprocess
import urllib.request
import tarfile
import zipfile
import shutil
from pathlib import Path


def get_system_info():
    """Get system information to determine the appropriate BLAST+ package."""
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


def check_if_installed():
    """Check if BLAST+ tools are already installed and accessible."""
    try:
        result = subprocess.run(["makeblastdb", "-version"], 
                              capture_output=True, text=True, check=True)
        print("✓ BLAST+ tools are already installed:")
        print(result.stdout.strip())
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def install_with_conda():
    """Try to install BLAST+ using conda if available."""
    try:
        # Check if conda is available
        subprocess.run(["conda", "--version"], capture_output=True, check=True)
        print("📦 Installing BLAST+ using conda...")
        
        # Install blast from bioconda
        result = subprocess.run([
            "conda", "install", "-c", "bioconda", "blast", "-y"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ Successfully installed BLAST+ using conda")
            return True
        else:
            print(f"⚠️  Conda installation failed: {result.stderr}")
            return False
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠️  Conda not available, trying manual installation...")
        return False


def install_with_homebrew():
    """Try to install BLAST+ using Homebrew on macOS."""
    try:
        # Check if homebrew is available
        subprocess.run(["brew", "--version"], capture_output=True, check=True)
        print("🍺 Installing BLAST+ using Homebrew...")
        
        result = subprocess.run([
            "brew", "install", "blast"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ Successfully installed BLAST+ using Homebrew")
            return True
        else:
            print(f"⚠️  Homebrew installation failed: {result.stderr}")
            return False
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠️  Homebrew not available, trying manual installation...")
        return False


def install_with_apt():
    """Try to install BLAST+ using apt on Ubuntu/Debian."""
    try:
        print("📦 Installing BLAST+ using apt...")
        
        # Update package list
        subprocess.run(["sudo", "apt", "update"], check=True)
        
        # Install ncbi-blast+
        result = subprocess.run([
            "sudo", "apt", "install", "-y", "ncbi-blast+"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ Successfully installed BLAST+ using apt")
            return True
        else:
            print(f"⚠️  Apt installation failed: {result.stderr}")
            return False
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠️  Apt not available, trying manual installation...")
        return False


def manual_install():
    """Manually download and install BLAST+ from NCBI."""
    print("📥 Downloading BLAST+ directly from NCBI...")
    
    try:
        system_type = get_system_info()
        
        # BLAST+ version - update this as needed
        blast_version = "2.15.0"
        
        # Determine download URL and file extension
        if system_type.startswith("win"):
            filename = f"ncbi-blast-{blast_version}+-{system_type}.exe"
            url = f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{blast_version}/{filename}"
        else:
            filename = f"ncbi-blast-{blast_version}+-{system_type}.tar.gz"
            url = f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{blast_version}/{filename}"
        
        print(f"Downloading {filename}...")
        
        # Create temporary directory
        temp_dir = Path("temp_blast_install")
        temp_dir.mkdir(exist_ok=True)
        
        # Download file
        file_path = temp_dir / filename
        urllib.request.urlretrieve(url, file_path)
        
        # Extract and install
        if system_type.startswith("win"):
            print("⚠️  Windows installer downloaded. Please run the .exe file manually.")
            print(f"File location: {file_path.absolute()}")
            return False
        else:
            print("📂 Extracting BLAST+ tools...")
            with tarfile.open(file_path, 'r:gz') as tar:
                tar.extractall(temp_dir)
            
            # Find the extracted directory
            extracted_dir = None
            for item in temp_dir.iterdir():
                if item.is_dir() and item.name.startswith("ncbi-blast"):
                    extracted_dir = item
                    break
            
            if not extracted_dir:
                raise RuntimeError("Could not find extracted BLAST+ directory")
            
            # Install to /usr/local/bin (or user's local bin)
            bin_dir = Path.home() / ".local" / "bin"
            bin_dir.mkdir(parents=True, exist_ok=True)
            
            blast_bin_dir = extracted_dir / "bin"
            for tool in blast_bin_dir.glob("*"):
                if tool.is_file():
                    shutil.copy2(tool, bin_dir / tool.name)
                    os.chmod(bin_dir / tool.name, 0o755)
            
            print(f"✓ BLAST+ tools installed to {bin_dir}")
            print(f"⚠️  Make sure {bin_dir} is in your PATH environment variable")
            
            # Add to PATH suggestion
            shell_rc = Path.home() / ".bashrc"
            if not shell_rc.exists():
                shell_rc = Path.home() / ".zshrc"
            
            print(f"\nTo add to PATH, add this line to your {shell_rc.name}:")
            print(f'export PATH="$HOME/.local/bin:$PATH"')
        
        # Cleanup
        shutil.rmtree(temp_dir)
        return True
        
    except Exception as e:
        print(f"❌ Manual installation failed: {e}")
        return False


def main():
    """Main installation function."""
    print("🧬 NCBI BLAST+ Installation Script")
    print("=" * 40)
    
    # Check if already installed
    if check_if_installed():
        return
    
    print("🔍 BLAST+ tools not found. Installing...")
    
    system = platform.system()
    
    # Try different installation methods based on the system
    success = False
    
    # First try conda (works on all systems)
    if install_with_conda():
        success = True
    
    # System-specific package managers
    elif system == "Darwin":  # macOS
        if install_with_homebrew():
            success = True
    elif system == "Linux":
        if install_with_apt():
            success = True
    
    # Fallback to manual installation
    if not success:
        print("📦 Trying manual installation...")
        success = manual_install()
    
    # Final check
    if success:
        print("\n🔄 Checking installation...")
        if check_if_installed():
            print("🎉 BLAST+ tools successfully installed!")
        else:
            print("⚠️  Installation completed but tools not found in PATH.")
            print("You may need to restart your terminal or update your PATH.")
    else:
        print("❌ Installation failed. Please install BLAST+ manually from:")
        print("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download")


if __name__ == "__main__":
    main()