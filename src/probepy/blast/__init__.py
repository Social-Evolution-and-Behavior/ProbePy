import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'install',
        'utils',
    },
    submod_attrs={
        'install': [
            'BLAST_VERSION',
            'check_blast_tools',
            'check_if_installed',
            'ensure_blast_tools',
            'get_system_info',
            'install_with_apt',
            'install_with_homebrew',
            'logger',
            'manual_blast_install',
        ],
        'utils': [
            'logger',
            'run_blastn',
            'run_makeblastdb',
        ],
    },
)

__all__ = ['BLAST_VERSION', 'check_blast_tools', 'check_if_installed',
           'ensure_blast_tools', 'get_system_info', 'install',
           'install_with_apt', 'install_with_homebrew', 'logger',
           'manual_blast_install', 'run_blastn', 'run_makeblastdb', 'utils']
