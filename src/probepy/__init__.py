"""
HCR-FISH Probe Design Toolkit

A comprehensive Python package for designing HCR-FISH probes compatible with
the HCR v3.0 amplifier system.
"""

# Import logging and set up package logger
import logging

# Create package logger
logger = logging.getLogger(__name__)

# Add a NullHandler by default so users can configure logging as they wish
logger.addHandler(logging.NullHandler())

# Import submodules for direct access
import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'blast',
        'hcr',
        'transcriptomics',
    },
    submod_attrs={
        'blast': [
            'BLAST_VERSION',
            'check_blast_tools',
            'check_if_installed',
            'ensure_blast_tools',
            'get_system_info',
            'install',
            'install_with_apt',
            'install_with_homebrew',
            'logger',
            'manual_blast_install',
            'run_blastn',
            'run_makeblastdb',
            'utils',
        ],
        'hcr': [
            'assign_target',
            'blast_gene',
            'check_probe_availability',
            'create_blast_databases',
            'design_hcr_probes',
            'download_with_rsync',
            'export_mrna_to_fasta',
            'get_amplifier',
            'get_probe_binding_regions_plot',
            'get_probes_IDT',
            'load_genome',
            'logger',
            'prep',
            'reverse_complement',
            'reverse_string',
            'set_logging_level',
            'utils',
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
            'logger',
            'main',
            'update',
            'update_transcriptome_object',
            'utils',
        ],
    },
)

__all__ = ['BLAST_VERSION', 'CDS', 'Exon', 'Gene', 'Intron', 'Transcript',
           'Transcriptome', 'UTR', 'assign_target', 'blast', 'blast_gene',
           'check_blast_tools', 'check_exons_contain_all_features',
           'check_if_installed', 'check_probe_availability', 'classes',
           'create_blast_databases', 'design_hcr_probes',
           'download_with_rsync', 'ensure_blast_tools', 'export_mrna_to_fasta',
           'generate_transcriptome_object', 'get_amplifier',
           'get_probe_binding_regions_plot', 'get_probes_IDT', 'get_sequence',
           'get_system_info', 'gtf_to_dataframe', 'hcr', 'install',
           'install_with_apt', 'install_with_homebrew', 'load_genome',
           'load_transcriptome_object', 'logger', 'main',
           'manual_blast_install', 'prep', 'reverse_complement',
           'reverse_string', 'run_blastn', 'run_makeblastdb',
           'set_logging_level', 'transcriptomics', 'update',
           'update_transcriptome_object', 'utils']
