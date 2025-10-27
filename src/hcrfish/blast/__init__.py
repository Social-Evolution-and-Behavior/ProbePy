import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'config',
        'install',
        'utils',
    },
    submod_attrs={
        'config': [
            'BLAST_VERSION',
        ],
        'install': [
            'check_blast_tools',
            'check_if_installed',
            'ensure_blast_tools',
            'get_system_info',
            'install_blast_tools',
            'install_with_apt',
            'install_with_homebrew',
            'manual_blast_install',
        ],
        'utils': [
            'run_blastn',
            'run_makeblastdb',
        ],
    },
)

__all__ = ['BLAST_VERSION', 'check_blast_tools', 'check_if_installed',
           'config', 'ensure_blast_tools', 'get_system_info', 'install',
           'install_blast_tools', 'install_with_apt', 'install_with_homebrew',
           'manual_blast_install', 'run_blastn', 'run_makeblastdb', 'utils']
