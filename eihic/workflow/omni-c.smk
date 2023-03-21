#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Snakemake for omni-c Hi-C libraries QC
"""
import os
import sys
import json
import csv

from eihic.scripts.hpc_config import HpcConfig
HPC_CONFIG = HpcConfig(config["hpc_config"])


# authorship
__author__ = "Mariano Ariel Olivera Fedi"
__maintainer__ = "Mariano Ariel Olivera Fedi"
__email__ = "mariano.olivera.fedi@hotmail.com"

# Request min version of snakemake
from snakemake.utils import min_version
min_version("7.10")

# Declare variables
jira_id = config["jira"]["jira_id"] # command line
#NOTIFY = not config["notify"] # command line
#if not Path(config["jira"]["password_file"]).is_file() or not jira_id:
#    NOTIFY = False

def list_basenames(paths):
    if paths is str:
        return os.path.basename(paths)
    
    basenames = []
    for path in paths:
        path = os.path.basename(path)
        basenames.append(path)
    return basenames


# Save the data to a list: R1 forward reads, R2 reverse reads, and sample name/ organism name
R1 = list_basenames(config["input_samples"]["R1"])
R2 = list_basenames(config["input_samples"]["R2"])
REFERENCE = list_basenames(config["input_samples"]["reference"])
ORGANISM = config["input_samples"]["organism"]
OUTPUT = config["output"]
logs = config["logs"]

if config['bwa-mem2'] == 'True':
    bwa = 'bwa-mem2'
    threading = "-p {threads}"
elif config['bwa-mem2'] == 'False':
    bwa = 'bwa',
    threading = ''
else:
    sys.exit("There is an error related to the config file.\n"
        + "Check your configuration file is not corrupted.\n"
        + "If you are not using eihic configure you should."
    )

# If one sample use string if more than one sample join the list.
# This is needed for one of the rules.
if len(R1) == 1 or len(R2) == 1:
    R1 = R1[0]
    R2 = R2[0]
else:
    R1 = " ".join(R1)
    R2 = " ".join(R2)


shell.prefix("set -eo pipefail; ")

rule all:
    input: 
        f"{OUTPUT}/results/{ORGANISM}_hi-c_library_complexity.txt", 
        f"{OUTPUT}/results/{ORGANISM}_stats_library.txt", 
        f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam"

rule library_complexity:
    input: f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam"
    output: f"{OUTPUT}/results/{ORGANISM}_hi-c_library_complexity.txt"
    log:
        os.path.join(logs, "library_complexity.log")
    params:
        source = config["source"]["omni-c"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u"
        + " && source {params.source}"
        + " && preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8" 
        + " -seg_len 1000000000 -output {output} {input}"
        + " ) > {log} 2>&1"

rule get_stats:
    input: f"{OUTPUT}/workflow/pairtools/stats.txt"
    output: f"{OUTPUT}/results/{ORGANISM}_stats_library.txt"
    log:
        os.path.join(logs, "get_stats.log")
    params:
        source = config["source"]["omni-c"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && get_qc.py -p {input} > {output}"
        + " ) > {log} 2>&1"


rule sort_bam:
    input: f"{OUTPUT}/workflow/pairtools/{ORGANISM}.bam"
    output: 
        sort = f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam", 
        index = f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam.bai"
    log:
        os.path.join(logs, "sort_bam.log")
    threads:
        int(HPC_CONFIG.get_cores("sort_bam"))
    params:
        source = config["source"]["omni-c"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u"
        + " && source {params.source}"
        + " && samtools sort -m 4G -@ {threads} -o {output.sort} {input}"
        + " && samtools index {output.sort}"
        + " ) > {log} 2>&1"

rule unsorted_bam:
    input: f"{OUTPUT}/workflow/pairtools/dedup.pairbam"
    output: 
        maps = f"{OUTPUT}/workflow/pairtools/mapped.pairs", 
        bam = temp(f"{OUTPUT}/workflow/pairtools/{ORGANISM}.bam")
    log:
        os.path.join(logs, "unsorted_bam.log")
    threads:
        int(HPC_CONFIG.get_cores("unsorted_bam"))/2
    params:
        source = config["source"]["omni-c"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u"
        + " && source {params.source}" 
        + " && pairtools split --nproc-in {threads} --nproc-out {threads}"
        + " --output-pairs {output.maps} --output-sam {output.bam} {input}"
        + " ) >  {log} 2>&1"

rule pairtools_dedup:
    input: f"{OUTPUT}/workflow/pairtools/sorted.pairbam"
    output: 
        out = temp(f"{OUTPUT}/workflow/pairtools/dedup.pairbam"), 
        stats =  f"{OUTPUT}/workflow/pairtools/stats.txt"
    log:
        os.path.join(logs, "pairtools_dedup.log")
    params:
        source = config["source"]["omni-c"]
    threads:
        int(HPC_CONFIG.get_cores("pairtools_dedup"))/2
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " &&  pairtools dedup --nproc-in {threads} --nproc-out {threads}"
        + " --mark-dups --output-stats {output.stats} --output {output.out} {input}"
        + " ) > {log} 2>&1"
    

rule pairtools_sort:
    input:f"{OUTPUT}/workflow/pairtools/parsed_pairbam"
    output:f"{OUTPUT}/workflow/pairtools/sorted.pairbam"
    params:
        source = config["source"]["omni-c"]
    threads:
        int(HPC_CONFIG.get_cores("pairtools_sort"))
    log:
        os.path.join(logs, "pairtools_sort.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + f" &&  pairtools sort --nproc {{threads}} --tmpdir={OUTPUT}/tmp {{input}} > {{output}}"
        + " ) > {log} 2>&1"
    

rule unique_unique:
    input: 
        align = f"{OUTPUT}/workflow/bwa/{ORGANISM}_mapped_reads.sort.bam",
        reference = f"{OUTPUT}/reference/genome/{ORGANISM}.genome"
    output: 
        unique = temp(f"{OUTPUT}/workflow/pairtools/parsed_pairbam"),
        stats = f"{OUTPUT}/workflow/pairtools/parsed_pairbam.stats"
    threads:
        HPC_CONFIG.get_cores("unique_unique")
    log:
        os.path.join(logs, "unique_unique.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    params: source = config["source"]["omni-c"]
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " &&  pairtools parse --min-mapq 40 --walks-policy 5unique" 
        + " --max-inter-align-gap 30 -o {output.unique}"
        + " --output-stats {output.stats}"
        + " --nproc-in {threads} --chroms-path {input.reference} {input.align}"
        + " ) > {log} 2>&1"


rule mapping_to_reference:
    input:
        R1=f"{OUTPUT}/reads/{R1}",
        R2=f"{OUTPUT}/reads/{R2}",
        bwa_inx = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta.amb",
        reference= f"{OUTPUT}/reference/genome/{ORGANISM}.fasta",
        index= f"{OUTPUT}/reference/genome/{ORGANISM}.fasta.fai"
    output: 
        f"{OUTPUT}/workflow/bwa/{ORGANISM}_mapped_reads.sort.bam"
    params: 
        source = config["source"]["omni-c"],
        bwa = bwa,
    threads:
        int(HPC_CONFIG.get_cores("mapping_to_reference"))
    log:
        os.path.join(logs, "mapping_to_reference.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " && {params.bwa} mem -5SP -T0 -t {threads} {input.reference}" 
        + " <(zcat -f {input.R1}) <(zcat -f {input.R2}) |"
        + " samtools view -@ {threads} -bS - |"
        + " samtools sort -m 4G -@ {threads} > {output}"
        + " ) > {log} 2>&1"

rule build_index:
    input: 
        reference = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta", 
        index = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta.fai", 
        genome = f"{OUTPUT}/reference/genome/{ORGANISM}.genome"
    output: 
        f"{OUTPUT}/reference/genome/{ORGANISM}.fasta.amb"
    params:
        source = config["source"]["omni-c"],
        bwa = bwa,
        threading = f"{threading}",
    threads:
        int(HPC_CONFIG.get_cores("build_index"))
    log:
        os.path.join(logs, "build_index.log")
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " && {params.bwa} index {params.threading} \ "
        + " {input.reference} || true ) > {log} 2>&1"

rule build_genome:
    input: 
        f"{OUTPUT}/reference/genome/{ORGANISM}.fasta.fai"
    output: 
        f"{OUTPUT}/reference/genome/{ORGANISM}.genome",
    log:
        os.path.join(logs, "build_genome.log")
    params:
        source = config["source"]["omni-c"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u"
        + " && source {params.source}" 
        + " && cut -f1,2 {input} > {output}"
        + " ) > {log} 2>&1"
    
rule index_reference:
    input: 
        f"{OUTPUT}/reference/genome/{ORGANISM}.fasta"      
    output: 
        f"{OUTPUT}/reference/genome/{ORGANISM}.fasta.fai"   
    log:
        os.path.join(logs, "index_reference.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    params:
        source = config["source"]["omni-c"]
    shell:
        "(set +u"
        + " && source {params.source}" 
        + " && samtools faidx {input} "
        + ") > {log} 2>&1"