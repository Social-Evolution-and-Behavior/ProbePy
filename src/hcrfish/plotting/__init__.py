import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'utils',
    },
    submod_attrs={
        'utils': [
            'black_plotting',
            'white_plotting',
        ],
    },
)

__all__ = ['black_plotting', 'utils', 'white_plotting']
