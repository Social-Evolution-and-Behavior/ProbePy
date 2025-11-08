import lazy_loader


__getattr__, __dir__, __all__ = lazy_loader.attach(
    __name__,
    submodules={
        'classes',
        'main',
        'update',
        'utils',
    },
    submod_attrs={
        'classes': [
            'CDS',
            'Exon',
            'Gene',
            'Intron',
            'Transcript',
            'Transcriptome',
            'UTR',
            'logger',
        ],
        'main': [
            'generate_transcriptome_object',
            'logger',
        ],
        'update': [
            'check_exons_contain_all_features',
            'load_transcriptome_object',
            'logger',
            'update_transcriptome_object',
        ],
        'utils': [
            'get_sequence',
            'gtf_to_dataframe',
        ],
    },
)

__all__ = ['CDS', 'Exon', 'Gene', 'Intron', 'Transcript', 'Transcriptome',
           'UTR', 'check_exons_contain_all_features', 'classes',
           'generate_transcriptome_object', 'get_sequence', 'gtf_to_dataframe',
           'load_transcriptome_object', 'logger', 'main', 'update',
           'update_transcriptome_object', 'utils']
