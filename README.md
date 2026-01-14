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

### Install Poetry 

Please follow the instructions [here](https://python-poetry.org/docs/) to install Poetry. 

### Install NCBI BLAST+ tools

Install [BLAST Command Line Tools 2.15.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/) from here. 
If running on a Mac, just use the installer `ncbi-blast-2.15.0+.dmg` 

**IMPORTANT**: Do not leave any spaces " " in the path to your project directory. It will break the BLAST toolkit. 


### Install ProbePy

Clone: 
```bash
git clone https://github.com/Social-Evolution-and-Behavior/ProbePy.git
cd ProbePy 
```

Optionally install the environment in your home directory: 

```bash
poetry config virtualenvs.in-project false
poetry config virtualenvs.path ~/.poetryvirtualenvs
```

Install: 
```bash
poetry lock
poetry install --with dev,notebook,docs -v
```

**Requirements:** Python 3.11+, Poetry, NCBI BLAST+ command line tools



## Quick Start

### 1. Basic Probe Design Workflow

```python
import probepy

# Check BLAST tools availability
probepy.check_blast_tools()

# Load transcriptome object (pre-built or create new)
transcriptome = probepy.load_transcriptome_object("species_transcriptome")

# Get gene of interest
gene = transcriptome.get_gene("Or9a")
transcript = gene.get_transcript_longest_cds()

# Design HCR probes
sequence = transcript.mrna_sequence
probes, regions, positions = probepy.design_hcr_probes(sequence, "B3")

print(f"Designed {len(probes)} probe pairs for {gene.name}")
```

### 2. Complete Pipeline Example

Visit the [docs](https://github.com/Social-Evolution-and-Behavior/ProbePy/tree/main/docs) for a complete walkthrough including:
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

## Output Formats

### IDT Order Sheets
- Excel format compatible with IDT oPools bulk ordering 
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


### Development 

After making changes please run 
```
poetry run mkinit --lazy_loader src/probepy --recursive --inplace 
```