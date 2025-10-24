"""
HCR-FISH Probe Design Toolkit

A comprehensive Python package for designing HCR-FISH probes compatible with
the HCR v3.0 amplifier system.
"""

# Import version information
from ._version import __version__

# Import submodules for direct access
from . import blast
from . import hcr  
from . import transcriptomics

__all__ = ['blast', 'hcr', 'transcriptomics', '__version__']
