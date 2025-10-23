#!/usr/bin/env python3
"""
Development utility script for the hcrfish package.

This script provides common development tasks like running tests,
formatting code, and checking code quality.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_tests(coverage=True, verbose=False):
    """Run the test suite."""
    print("🧪 Running tests...")
    
    cmd = ["python", "-m", "pytest"]
    
    if coverage:
        cmd.extend(["--cov=hcrfish", "--cov-report=term-missing"])
    
    if verbose:
        cmd.append("-v")
    
    result = subprocess.run(cmd)
    return result.returncode == 0


def format_code():
    """Format code with black and isort."""
    print("🎨 Formatting code...")
    
    # Run black
    print("  Running black...")
    black_result = subprocess.run(["python", "-m", "black", "src", "tests"])
    
    # Run isort
    print("  Running isort...")
    isort_result = subprocess.run(["python", "-m", "isort", "src", "tests"])
    
    return black_result.returncode == 0 and isort_result.returncode == 0


def check_code_quality():
    """Check code quality with flake8 and mypy."""
    print("🔍 Checking code quality...")
    
    # Run flake8
    print("  Running flake8...")
    flake8_result = subprocess.run(["python", "-m", "flake8", "src", "tests"])
    
    # Run mypy (may fail if not all dependencies have types)
    print("  Running mypy...")
    mypy_result = subprocess.run(["python", "-m", "mypy", "src"], 
                               capture_output=True, text=True)
    
    if mypy_result.returncode != 0:
        print("  ⚠️  MyPy warnings (some dependencies may lack type annotations):")
        print(mypy_result.stdout)
    
    return flake8_result.returncode == 0


def check_blast_tools():
    """Check if BLAST tools are available."""
    print("🔬 Checking BLAST tools...")
    
    try:
        import hcrfish.blast
        status = hcrfish.blast.check_blast_tools()
        
        all_available = all(tool['available'] for tool in status.values())
        
        if all_available:
            print("  ✅ All BLAST tools are available")
            return True
        else:
            print("  ❌ Some BLAST tools are missing")
            for tool, info in status.items():
                if not info['available']:
                    print(f"    - {tool}: Not found")
            return False
            
    except ImportError as e:
        print(f"  ❌ Could not import hcrfish.blast: {e}")
        return False


def install_dev_dependencies():
    """Install development dependencies."""
    print("📦 Installing development dependencies...")
    
    result = subprocess.run([
        "python", "-m", "pip", "install", 
        "black", "isort", "flake8", "mypy", "pytest", "pytest-cov", "pytest-mock"
    ])
    
    return result.returncode == 0


def main():
    """Main CLI interface."""
    parser = argparse.ArgumentParser(description="Development utilities for hcrfish")
    
    parser.add_argument("--test", action="store_true", 
                       help="Run the test suite")
    parser.add_argument("--no-coverage", action="store_true",
                       help="Skip coverage reporting when running tests")
    parser.add_argument("--format", action="store_true",
                       help="Format code with black and isort")
    parser.add_argument("--check", action="store_true",
                       help="Check code quality with flake8 and mypy")
    parser.add_argument("--blast", action="store_true",
                       help="Check BLAST tools availability")
    parser.add_argument("--install-dev", action="store_true",
                       help="Install development dependencies")
    parser.add_argument("--all", action="store_true",
                       help="Run all checks (format, test, quality)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    # If no specific action, show help
    if not any([args.test, args.format, args.check, args.blast, 
                args.install_dev, args.all]):
        parser.print_help()
        return 0
    
    success = True
    
    if args.install_dev:
        success &= install_dev_dependencies()
    
    if args.all or args.format:
        success &= format_code()
    
    if args.all or args.test:
        success &= run_tests(coverage=not args.no_coverage, verbose=args.verbose)
    
    if args.all or args.check:
        success &= check_code_quality()
    
    if args.blast:
        success &= check_blast_tools()
    
    if success:
        print("\n✅ All checks passed!")
        return 0
    else:
        print("\n❌ Some checks failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())