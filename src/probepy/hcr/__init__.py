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
            'logger',
        ],
        'utils': [
            'assign_target',
            'blast_gene',
            'check_probe_availability',
            'design_hcr_probes',
            'get_amplifier',
            'get_probe_binding_regions_plot',
            'get_probes_IDT',
            'load_genome',
            'logger',
            'reverse_complement',
            'reverse_string',
            'set_logging_level',
        ],
    },
)

__all__ = ['assign_target', 'blast_gene', 'check_probe_availability',
           'create_blast_databases', 'design_hcr_probes',
           'download_with_rsync', 'export_mrna_to_fasta', 'get_amplifier',
           'get_probe_binding_regions_plot', 'get_probes_IDT', 'load_genome',
           'logger', 'prep', 'reverse_complement', 'reverse_string',
           'set_logging_level', 'utils']
