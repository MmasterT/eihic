# Eihic

This tool is an in-house snakemake pipeline from the EI Core-Bioinformatics group for the quality control and pre-processinf of Hi-C data for genome assembly. This tool currently suports the Arima Genomics and Dovetail-Genomics Omni-C protocols.

## Installation

```bash
git clone https://github.com/EI-CoreBioinformatics/eihic.git
cd eihic/eihic
pip install .
```

## Quick start

## Documentation

### Omni-c

[https://omni-c.readthedocs.io/en/latest/](https://omni-c.readthedocs.io/en/latest/)

https://pairtools.readthedocs.io/en/latest/

https://bio-bwa.sourceforge.net/bwa.shtml

### Arima

https://github.com/ArimaGenomics/CHiC

HiCUP is the software used for the QC from

[https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html](https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html)

[https://bioconductor.org/packages/release/bioc/html/Chicago.html](https://bioconductor.org/packages/release/bioc/html/Chicago.html)

The installation of the Arima software is sourced from a fork of the original repository because of some subtle changes to the installation steps were necessary for the HPC installation.

https://github.com/MmasterT/CHiC
