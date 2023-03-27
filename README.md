# EIHiC

This tool is an in-house snakemake pipeline from the EI Core-Bioinformatics group for the quality control and pre-processinf of Hi-C data for genome assembly. This tool currently suports the Arima Genomics and Dovetail-Genomics Omni-C protocols.

## Installation

```bash
#Snakemake needs to be at least version 7.10.
#This software was tested with 7.12.1 and it is recommended to use the same due major and inconsistent changes in snakemake.

conda create snakemake python=3.8 snakemake=7.12.1
conda activate snakemake #If it is your first time go to the conda documentation at the end of the repository. 

git clone https://${PERSONAL_TOKEN}@github.com/EI-CoreBioinformatics/eihic.git

input USERNAME
input PERSONAL_TOKEN

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

First source the current version of the tool:

`source eihic-0.2.0`

If you don't have it add /ei/software/cb/bin to your PATH variable:

`export PATH=${PATH}:/ei/software/cb/bin`

Then check the wrapper of the two main companents of eihic:

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

The first step of the pipeline is to create a config file to be used by snakemake in the HPC. All the arguments are described as optional but --sample_csv -s  must be have an asociated Path a sample file as showed in the reference path in the eihic configure --help command.

```bash
--% eihic configure --help
usage: EI HI-C configure [-h] [-s SAMPLES_CSV] [-c] [--jira JIRA] [-o OUTPUT] [-f] [-bm2]

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLES_CSV, --samples_csv SAMPLES_CSV
                        Provide sample information in tab-separated format. Please refer to the sample
                        file: /ei/software/cb/eihic/0.2.0/x86_64/lib/python3.9/site-
                        packages/eihic/etc/run_config.yaml for more information above the csv format. A
                        template is provided here
                        /ei/software/cb/eihic/0.2.0/x86_64/lib/python3.9/site-
                        packages/eihic/etc/samples.csv. (default: None)
  -c, --curation        Use this flag if you have your final scaffold and want to run the curation
                        pipeline. Bare in mind that the sample_csv file requires an extra field with
                        the long reads file paths as a fifth line. (default: False)
  --jira JIRA           Provide JIRA id for posting job summary. E.g., PPBFX-611 (default: None)
  -o OUTPUT, --output OUTPUT
                        Provide output directory (default: /ei/.project-
                        scratch/3/35652c26-5607-45bf-a6ec-2fd255ce2730/CB-
                        GENANNO-532_ERGA_Salvelinus_alpinus/output)
  -f, --force-reconfiguration
                        Force reconfiguration (default: False)
  -bm2, --bwa_mem2      Use bwa-mem2 insted of bwa mem. This option use a lot more RAM, use with
                        precaution. (default: False
```

Once you created the output_dir/run_config.yaml file from the eihic wrapper you need to run the eihic run subcommand. This tool supports omni-c and arima two-enzymes library prep protocols fot Hi-C, you must choose one of them (default is omni-c).

```bash
--% eihic run --help
usage: EI HI-C run [-h] [--library LIBRARY] [-c] [--hpc_config HPC_CONFIG] [--jobs JOBS]
                   [--latency_wait LATENCY_WAIT] [--no_posting_off] [-v] [-np]
                   run_config

positional arguments:
  run_config            Provide run configuration YAML. Run 'eihic configure -h' to generate the run
                        configuration YAML file. (Description template file is here:
                        /ei/software/cb/eihic/0.2.0/x86_64/lib/python3.9/site-
                        packages/eihic/etc/run_config.yaml)

optional arguments:
  -h, --help            show this help message and exit
  --library LIBRARY     This pipeline supports the following library protocols for hi-c data: arima (2
                        enzymes protocol), and omni-c. Provide the name (arima or omni-c) as the second
                        positional argument after the sample file (default: omni-c)
  -c, --curation        This flag run the steps generate all the required Hi-C contact matrices and
                        tracks for its manual curation step. (default: False)
  --hpc_config HPC_CONFIG
                        Provide HPC configuration YAML (default:
                        /ei/software/cb/eihic/0.2.0/x86_64/lib/python3.9/site-
                        packages/eihic/etc/hpc_config.json)
  --jobs JOBS, -j JOBS  Use at most N CPU cluster/cloud jobs in parallel (default: 100)
  --latency_wait LATENCY_WAIT
                        Wait given seconds if an output file of a job is not present after the job
                        finished (default: 120)
  --no_posting_off      Use this flag if you want to post comments to JIRA tickets (default: True)
  -v, --verbose         Verbose mode for debugging (default: False)
  -np, --dry_run        Dry run (default: False
```

## Relevant outputs of each pipeline:



arima two enzymes library:

- READ_NAME.hicup.sam.HiCUP_summary_report_UUID_TIMESTAMP.html


omni-c library:

- OUTPUT_dir/results/ORGANISM_NAME_hi-c_library_complexity.txt",

- OUTPUT_dir/results/ORGANISM_NAME_stats_library.txt",

- OUTPUT_dir/workflow/samtools/ORGANISM_NAME.sorted.bam

- OUTPUT_dir/bwa/samtools/ORGANISM_NAME.sorted.bam


omni-c library curation mode:

- OUTPUT_dir/workflow/samtools/ORGANISM_NAME.sorted.bam
- OUTPUT_dir/workflow/bwa/ORGANISM_NAME_mapped_reads.sort.bam
- OUTPUT_dir/workflow/pretext/ORGANISM_NAME_unique_mapping.pretext
- OUTPUT_dir/workflow/pretext/ORGANISM_NAME_multi_mapping.pretext
- OUTPUT_dir/workflow/cooler/unique_1kb.mcool
- OUTPUT_dir/workflow/cooler/all_1kb.mcool
- OUTPUT_dir/workflow/tracks/gaps_ORGANISM_NAME.bedgraph
- OUTPUT_dir/workflow/tracks/telomeres_ORGANISM_NAME.bedgraph
- OUTPUT_dir/workflow/tracks/coverage_ORGANISM_NAME.bedgraph

## Documentation

### Sample_csv format



This file is composed of four lines:

1. all R1 reads
2. all R2 reads
3. Path to the reference assembly.
4. organism name (will be used for file naming)
5. (if running curation mode -c / --curation) list of hifi_reads in fasta/fasta.gz



Example:

sample_1_R1.fastq,sample_2_R1.fastq,(...), sample_n_R1.fastq

sample_1_R2.fastq,sample_2_R2.fastq,(...), sample_n_R2.fastq

reference/genome/you_reference_genome.fasta #path to reference

name_of_organism (internal naming usage)

hifi_1.fasta,hifi_2.fasta,(...),hifi_n.fasta (optional)

### Conda

https://docs.conda.io/projects/conda/en/stable/

Package managment and environment mangament for python code.

### SNAKEMAKE

Workflow managment tool based in Python.

https://snakemake.readthedocs.io/en/v7.12.0/

https://github.com/snakemake/snakemake

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

### Pairix

https://github.com/4dn-dcic/pairix

Matrix indexing used to obtain .mcool contac maps to use in high-glass

### Cooler

https://github.com/open2c/cooler

https://cooler.readthedocs.io/en/latest/

multiple resolution matrix used for the manual curation of genomes and TAD analysis.

### PretextMap

https://github.com/wtsi-hpag/PretextMap

Contact maps visualization tool for Hi-C. Used to move contigs/scaffolds to curate the genome.

### TIDK

Telomere identification toolkit from DToL.

https://github.com/tolkit/telomeric-identifier

If you don't know or there is no information for your telomere sequence

### Mosdepth

https://github.com/brentp/mosdepth

obtaining coverage data from the genome.

### Minimap2

https://github.com/lh3/minimap2

aligner of choice for the PacBio HiFi reads. The output is in sam/bam format and uses the the preset map-hifi

## To Do's

- use bwa-mem 2 to make analysis fastaer for small genomes (too resource demanding for big genomes)
- add other relevant library preps for hi-c
- add telomere to the config
- add repeats trach from eirepeat pipeline
-
