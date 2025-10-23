import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'utils',
    },
    submod_attrs={
        'utils': [
            'design_hcr_probes',
            'get_amplifier',
            'reverse_complement',
            'reverse_string',
        ],
    },
)

__all__ = ['design_hcr_probes', 'get_amplifier', 'reverse_complement',
           'reverse_string', 'utils']
