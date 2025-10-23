"""
HCR-FISH Probe Design Toolkit

A comprehensive Python package for designing HCR-FISH probes compatible with
the HCR v3.0 amplifier system.
"""

import lazy_loader

# Import version information
from ._version import __version__, __author__, __email__, __description__

__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'blast',
        'hcr',
        'plotting',
        'transcriptomics',
    },
    submod_attrs={
        'blast': [
            'blast_utils',
            'check_blast_tools',
            'check_if_installed',
            'ensure_blast_tools',
            'get_system_info',
            'install_ncbi_tools',
            'install_with_apt',
            'install_with_conda',
            'install_with_homebrew',
            'main',
            'manual_install',
            'run_blastn',
            'run_makeblastdb',
        ],
        'hcr': [
            'design_hcr_probes',
            'get_amplifier',
            'reverse_complement',
            'reverse_string',
            'utils',
        ],
        'plotting': [
            'black_plotting',
            'utils',
            'white_plotting',
        ],
        'transcriptomics': [
            'CDS',
            'Exon',
            'Gene',
            'Intron',
            'Transcript',
            'Transcriptome',
            'UTR',
            'check_exons_contain_all_features',
            'classes',
            'generate_transcriptome_object',
            'get_sequence',
            'gtf_to_dataframe',
            'load_transcriptome_object',
            'main',
            'update',
            'update_transcriptome_object',
            'utils',
        ],
    },
)

__all__ = ['CDS', 'Exon', 'Gene', 'Intron', 'Transcript', 'Transcriptome',
           'UTR', 'black_plotting', 'blast', 'blast_utils',
           'check_blast_tools', 'check_exons_contain_all_features',
           'check_if_installed', 'classes', 'design_hcr_probes',
           'ensure_blast_tools', 'generate_transcriptome_object',
           'get_amplifier', 'get_sequence', 'get_system_info',
           'gtf_to_dataframe', 'hcr', 'install_ncbi_tools', 'install_with_apt',
           'install_with_conda', 'install_with_homebrew',
           'load_transcriptome_object', 'main', 'manual_install', 'plotting',
           'reverse_complement', 'reverse_string', 'run_blastn',
           'run_makeblastdb', 'transcriptomics', 'update',
           'update_transcriptome_object', 'utils', 'white_plotting']
