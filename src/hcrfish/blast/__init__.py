import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'blast_utils',
        'install_ncbi_tools',
    },
    submod_attrs={
        'blast_utils': [
            'check_blast_tools',
            'ensure_blast_tools',
            'run_blastn',
            'run_makeblastdb',
        ],
        'install_ncbi_tools': [
            'check_if_installed',
            'get_system_info',
            'install_with_apt',
            'install_with_conda',
            'install_with_homebrew',
            'main',
            'manual_install',
        ],
    },
)

__all__ = ['blast_utils', 'check_blast_tools', 'check_if_installed',
           'ensure_blast_tools', 'get_system_info', 'install_ncbi_tools',
           'install_with_apt', 'install_with_conda', 'install_with_homebrew',
           'main', 'manual_install', 'run_blastn', 'run_makeblastdb']
