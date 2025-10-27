import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'prep',
        'utils',
    },
    submod_attrs={
        'prep': [
            'create_blast_databases',
            'download_with_rsync',
            'export_mrna_to_fasta',
        ],
        'utils': [
            'blast_gene',
            'check_probe_availability',
            'design_hcr_probes',
            'get_probe_binding_regions_plot',
            'get_amplifier',
            'get_probes_IDT',
            'reverse_complement',
            'reverse_string',
        ],
    },
)

__all__ = ['blast_gene', 'check_probe_availability', 'create_blast_databases',
           'design_hcr_probes', 'download_with_rsync', 'export_mrna_to_fasta',
           'get_probe_binding_regions_plot', 'get_amplifier',
           'get_probes_IDT', 'prep', 'reverse_complement', 'reverse_string',
           'utils']
