"""
HCR-FISH Probe Design Toolkit

A comprehensive Python package for designing HCR-FISH probes compatible with
the HCR v3.0 amplifier system.
"""

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
            'utils',
        ],
        'hcr': [
            'blast_gene',
            'check_probe_availability',
            'create_blast_databases',
            'design_hcr_probes',
            'download_with_rsync',
            'export_mrna_to_fasta',
            'get_probe_binding_regions_plot',
            'get_amplifier',
            'get_probes_IDT',
            'prep',
            'reverse_complement',
            'reverse_string',
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
            'main',
            'update',
            'update_transcriptome_object',
            'utils',
        ],
    },
)

__all__ = ['CDS', 'Exon', 'Gene', 'Intron', 'Transcript', 'Transcriptome',
           'UTR', 'blast', 'blast_gene', 'check_blast_tools',
           'check_exons_contain_all_features', 'check_if_installed',
           'check_probe_availability', 'classes', 'create_blast_databases',
           'design_hcr_probes', 'download_with_rsync', 'ensure_blast_tools',
           'export_mrna_to_fasta', 'get_probe_binding_regions_plot',
           'generate_transcriptome_object', 'get_amplifier', 'get_probes_IDT',
           'get_sequence', 'get_system_info', 'gtf_to_dataframe', 'hcr',
           'install_ncbi_tools', 'install_with_apt', 'install_with_conda',
           'install_with_homebrew', 'load_transcriptome_object', 'main',
           'manual_install', 'prep', 'reverse_complement', 'reverse_string',
           'run_blastn', 'run_makeblastdb', 'transcriptomics', 'update',
           'update_transcriptome_object', 'utils']
