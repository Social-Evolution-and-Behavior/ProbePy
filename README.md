## ProbePy

> **Fluorescent *in situ* Hybridization Probe Design Tool**

A comprehensive Python package for designing custom fluorescent *in situ* hybridization (FISH) probes. 

## Features

- **Multi-species support**: Works with any organism with genome and transcriptome annotations
- **Compatible probe design**: Probes are compatible with Molecular Instruments HCR v3.0 B-system amplification 
- **Facilitates Multiplexing**: Probes can be ordered with B1-B5 amplifier recognition motifs 
- **BLAST-based specificity**: Ensures unique binding by blasting against the transcriptome 
- **Comprehensive visualization**: Generate publication-ready plots showing probe locations on gene structures
- **Flexible output formats**: Export probes in IDT-compatible Excel format for direct ordering

## Installation

```bash
git clone https://github.com/Social-Evolution-and-Behavior/ProbePy.git
cd ProbePy 
poetry install
```

**Requirements:** Python 3.11+, Poetry, NCBI BLAST+ command line tools

### Install NCBI BLAST+ tools

The package facilitates the installation of NCBI BLAST+ command line tools. 

```python
import hcrfish
hcrfish.install_blast_tools()
```

Alternatively, you can install it yourself: 
- **macOS**: `brew install blast`
- **Ubuntu/Debian**: `sudo apt-get install ncbi-blast+`
- **Windows**: Download from [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)

## Quick Start

### 1. Basic Probe Design Workflow

```python
import hcrfish

# Check BLAST tools availability
hcrfish.check_blast_tools()

# Load transcriptome object (pre-built or create new)
transcriptome = hcrfish.load_transcriptome_object("species_transcriptome")

# Get gene of interest
gene = transcriptome.get_gene("Or9a")
transcript = gene.get_transcript_longest_cds()

# Design HCR probes
sequence = transcript.mrna_sequence
probes, regions, positions = hcrfish.design_hcr_probes(sequence, "B3")

print(f"Designed {len(probes)} probe pairs for {gene.name}")
```

### 2. Complete Pipeline Example

Visit the [docs] for a complete walkthrough including:
- Downloading data with rsync
- Building BLAST databases 
- Transcriptome loading and gene selection
- BLAST-based specificity checking
- Probe design and filtering
- Visualization and export


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


## Output Formats

### IDT Order Sheets
- Excel format compatible with IDT bulk ordering
- Includes pool names and probe sequences
- Ready for direct upload to IDT website

### Probe Binding Regions
- Detailed information about where each probe binds
- Useful for troubleshooting and analysis

### Visualization
- Gene structure plots with probe locations

## Acknowledgments

- HCR v3.0 system developed by Molecular Instruments 
- Built with Python scientific computing stack (NumPy, Pandas, Matplotlib)
- BLAST integration via NCBI BLAST+ suite
- Visualization powered by pygenomeviz 


