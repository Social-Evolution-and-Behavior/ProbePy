# HCR-FISH Probe Design Tool

[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A comprehensive Python package for designing HCR-FISH (Hybridization Chain Reaction Fluorescent In Situ Hybridization) probes compatible with the HCR v3.0 amplifier system. This tool enables researchers to design highly specific, custom probes for any species with available genome and transcriptome data.

## Features

- **Automated probe design**: Design HCR v3.0 compatible probe pairs with optimal GC content and minimal off-target binding
- **Multi-species support**: Works with any organism with genome and transcriptome annotations
- **BLAST-based specificity verification**: Ensures probe specificity using NCBI BLAST+ tools
- **Comprehensive visualization**: Generate publication-ready plots showing probe locations on gene structures
- **Flexible output formats**: Export probes in IDT-compatible Excel format for direct ordering
- **Built-in quality control**: Automatic filtering for problematic sequences (repeats, extreme GC content)

## Installation

### Prerequisites

- Python 3.11 or higher
- Poetry (for dependency management)
- NCBI BLAST+ command line tools

### Quick Setup

1. **Clone the repository:**
```bash
git clone <repository-url>
cd hcr_probe_design_general
```

2. **Install dependencies with Poetry:**
```bash
poetry install
poetry shell
```

3. **Install NCBI BLAST+ tools:**
The package will automatically attempt to install BLAST+ tools when first used. For manual installation:
- **macOS**: `brew install blast`
- **Ubuntu/Debian**: `sudo apt-get install ncbi-blast+`
- **Windows**: Download from [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)

## Quick Start

### 1. Basic Probe Design Workflow

```python
import hcrfish

# Check BLAST tools availability
hcrfish.blast.check_blast_tools()

# Load transcriptome object (pre-built or create new)
transcriptome = hcrfish.load_transcriptome_object("species_transcriptome")

# Get gene of interest
gene = transcriptome.get_gene("Or9a")
transcript = gene.get_transcript_longest_bounds()

# Design HCR probes
sequence = transcript.mrna_sequence
probes, regions, positions = hcrfish.design_hcr_probes(sequence, "B3")

print(f"Designed {len(probes)} probe pairs for {gene.name}")
```

### 2. Complete Pipeline Example

See `design-probes.ipynb` for a complete walkthrough including:
- Transcriptome loading and gene selection
- BLAST-based specificity checking
- Probe design and filtering
- Visualization and export

## Core Modules

### `hcrfish.transcriptomics`
- **Classes**: `Transcriptome`, `Gene`, `Transcript`, `Exon`, `CDS`, `UTR`, `Intron`
- **Functions**: Build and manage transcriptome objects from GTF/GFF files

### `hcrfish.hcr`
- **`design_hcr_probes()`**: Core probe design algorithm
- **`get_amplifier()`**: HCR v3.0 amplifier sequences (B1-B5)
- **`reverse_complement()`**: DNA sequence utilities

### `hcrfish.blast`
- **`check_blast_tools()`**: Verify BLAST+ installation
- **`run_blastn()`**: Automated BLAST searches for specificity
- **`run_makeblastdb()`**: Database creation utilities

### `hcrfish.plotting`
- **`white_plotting()`** / **`black_plotting()`**: Publication-ready plot themes
- Integration with `pygenomeviz` for gene structure visualization

## HCR v3.0 Amplifier Support

The package supports all HCR v3.0 amplifiers:

| Amplifier | Upstream Initiator | Downstream Initiator |
|-----------|-------------------|---------------------|
| B1 | GAGGAGGGCAGCAAACGG | GAAGAGTCTTCCTTTACG |
| B2 | CCTCGTAAATCCTCATCA | ATCATCCAGTAAACCGCC |
| B3 | GTCCCTGCCTCTATATCT | CCACTCAACTTTAACCCG |
| B4 | CCTCAACCTACCTCCAAC | TCTCACCATATTCGCTTC |
| B5 | CTCACTCCCAATCTCTAT | CTACCCTACAAATCCAAT |

## Probe Design Algorithm

The probe design follows these key principles:

1. **Target sequence processing**: Uses mature mRNA (exons only, no introns)
2. **Probe pair generation**: Designs complementary 25-nt probe pairs with optimal spacing
3. **Quality filtering**: 
   - GC content between 25-75%
   - No homopolymer runs >4 bases
   - No gaps from unique region masking
4. **Specificity verification**: BLAST against transcriptome with/without introns
5. **Visualization**: Maps probe locations to genomic coordinates

## Example Workflows

### Design probes for a single gene
```python
# Load transcriptome
tr = hcrfish.load_transcriptome_object("mel_transcriptome")

# Select gene and transcript
gene = tr.get_gene("dsx")
transcript = gene.get_transcript_longest_bounds()

# Design 30 B3 probes
probes, regions, positions = hcrfish.design_hcr_probes(transcript.mrna_sequence, "B3")
selected_probes = probes[:30]  # Take first 30 probes
```

### Batch processing multiple genes
```python
genes = ["Or9a", "dsx", "fru"]
all_probes = {}

for gene_name in genes:
    gene = tr.get_gene(gene_name)
    transcript = gene.get_transcript_longest_bounds()
    probes, regions, positions = hcrfish.design_hcr_probes(transcript.mrna_sequence, "B3")
    all_probes[gene_name] = probes[:20]  # 20 probes per gene
```

## Output Formats

### IDT Order Sheets
- Excel format compatible with IDT bulk ordering
- Includes pool names and probe sequences
- Ready for direct upload to IDT website

### Probe Binding Regions
- Detailed information about where each probe binds
- Coordinates in both transcript and genomic space
- Useful for troubleshooting and analysis

### Visualization
- Gene structure plots with probe locations
- Publication-ready PNG exports
- Both light and dark themes available

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

### Development Setup

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Install development dependencies: `poetry install --with dev`
4. Run tests: `pytest`
5. Submit pull request

## Citation

If you use this tool in your research, please cite:

```
HCR-FISH Probe Design Tool
Giacomo Glotzer, The Rockefeller University
https://github.com/your-username/hcr_probe_design_general
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Contact: gglotzer@rockefeller.edu

## Acknowledgments

- HCR v3.0 system developed by Molecular Instruments
- Built with Python scientific computing stack (NumPy, Pandas, Matplotlib)
- BLAST integration via NCBI BLAST+ suite
- Visualization powered by pygenomeviz 


