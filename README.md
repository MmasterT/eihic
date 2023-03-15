# EIHiC

This tool is an in-house snakemake pipeline from the EI Core-Bioinformatics group for the quality control and pre-processinf of Hi-C data for genome assembly. This tool currently suports the Arima Genomics and Dovetail-Genomics Omni-C protocols.

## Installation

```bash
#Snakemake needs to be at least version 7.10.
#This software was tested with 7.12.1 and it is recommended to use the same due major and inconsistent changes in snakemake.

conda create snakemake python=3.8 snakemake=7.12.1
conda activate snakemake #If it is your first time go to the conda documentation at the end of the repository. 

git clone https://github.com/EI-CoreBioinformatics/eihic.git
cd eihic

#For ease of installation follow these steps

python -m pip install --upgrade pip #updates pip to latest version

PYTHONUSERBASE=/path/to/install/to #Sets a custom path of installation
python3 -m pip install --user -r requirements.txt .

echo $PATH #be sure the binaries are in your Path

#if not add to path.
PATH="${PYTHONUSERBASE}"/bin:"${PATH}" 

```

## Quick start

Let's check the wrapper of the two main companents of EIHiC:

```bash
Usage:
 eihic  --help

[~]--% eihic --help
usage: EI HI-C [-h] [-v] {configure,run} ...

EI HI-C Pipeline

positional arguments:
  {configure,run}
    configure      see `configure -h`
    run            see `run -h`

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show program's version number and exit
```

The first step of the pipeline is to create a config file to be used by snakemake in the HPC. All the arguments are described as optional but --sample_csv must be have an asociated Path ro a sample file as showed in the reference path in the eihic configure --help command.

```bash
[~]--% eihic configure --help
usage: EI HI-C configure [-h] [--samples_csv SAMPLES_CSV] [--jira JIRA] [-o OUTPUT] [-f] [-bm2]

optional arguments:
  -h, --help            show this help message and exit
  --samples_csv SAMPLES_CSV
                        Provide sample information in tab-separated format. Please refer to the sample file:
                        /ei/software/cb/eihic/0.1.0/x86_64/lib/python3.9/site-packages/eihic/etc/run_config.yaml, for more information above the csv format. A template is provided here
                        /ei/software/cb/eihic/0.1.0/x86_64/lib/python3.9/site-packages/eihic/etc/samples.csv.
                        (default: None)
  --jira JIRA           Provide JIRA id for posting job summary. E.g., PPBFX-611 (default: None)
  -o OUTPUT, --output OUTPUT
                        Provide output directory (default: /hpc-home/olivera/output)
  -f, --force-reconfiguration
                        Force reconfiguration (default: False)
  -bm2, --bwa_mem2      Use bwa-mem2 insted of bwa mem. This option use a lot more RAM, use with precaution.
                        (default: False)
```

Once you got the run_config.yaml file from the EIHiC wrapper you need to specify the run command. This tool supports omni-c and arima two-enzymes library prep protocols fot Hi-C, you must choose one of them (default is omni-c).

```bash
[~]--% eihic run --help
usage: EI HI-C run [-h] [--library LIBRARY] [--hpc_config HPC_CONFIG] [--jobs JOBS] [--latency_wait LATENCY_WAIT]
                   [--no_posting] [-v] [-np]
                   run_config

positional arguments:
  run_config            Provide run configuration YAML. Run 'eihic configure -h' to generate the run configuration
                        YAML file. (Description template file is here:
                        /ei/software/cb/eihic/0.1.0/x86_64/lib/python3.9/site-packages/eihic/etc/run_config.yaml)

optional arguments:
  -h, --help            show this help message and exit
  --library LIBRARY     This pipeline supports the following library protocols for hi-c data: arima (2 enzymes
                        protocol), and omni-c. Provide the name (arima or omni-c) as the second positional argument
                        after the sample file (default: omni-c).
  --hpc_config HPC_CONFIG
                        Provide HPC configuration YAML (default:
                        /ei/software/cb/eihic/0.1.0/x86_64/lib/python3.9/site-packages/eihic/etc/hpc_config.json)
  --jobs JOBS, -j JOBS  Use at most N CPU cluster/cloud jobs in parallel (default: 100)
  --latency_wait LATENCY_WAIT
                        Wait given seconds if an output file of a job is not present after the job finished
                        (default: 120)
  --no_posting_off          Use this flag if you are testing and do not want to post comments to JIRA tickets (default:
                        True)
  -v, --verbose         Verbose mode for debugging (default: False)
  -np, --dry_run        Dry run (default: False)
```

## Documentation

### Sample_csv format

This file is composed of four lines:

1. all R1 reads
2. all R2 reads
3. Path to the reference assembly.
4. organism name (will be used for file naming)
5. (if running curation mode) list of hifi_reads in fasta/fasta.gz

Example:

sample_1_R1.fastq,sample_2_R1.fastq,(...), sample_n_R1.fastq

sample_1_R2.fastq,sample_2_R2.fastq,(...), sample_n_R2.fastq

reference/genome/you_reference_genome.fasta #path to reference

name_of_organism (internal naming usage)

hifi_1.fasta,hifi_2.fasta,(...),hifi_n.fasta

### Conda

https://docs.conda.io/projects/conda/en/stable/

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

### Pairbix

### Cooler

### PretextMap

### TIDK

### Mosdepth

### YaHS

### Minimap2

## To Do's

- find how to delete the bug that does not let you use absoulth path for the files
- use bwa-mem 2 to make analysis fastaer for small genomes (too resource demanding for big genomes)
- add other relevant library preps for hi-c
