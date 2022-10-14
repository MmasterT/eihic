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
NOTIFY = not config["notify"] # command line
if not Path(config["jira"]["password_file"]).is_file() or not jira_id:
    NOTIFY = False

# Open sample_information.csv and parse the data
with open('reference/sample_information.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

# Save the data to a list: R1 forward reads, R2 reverse reads, and sample name/ organism name
R1= config["R1"]
R2= config["R2"]
REFERENCE= config["reference"]
ORGANISM= config["organism"]

# If one sample use string if more than one sample join the list
if len(data[0]) == 1 or len(data[1]) == 1:
    R1 = data[0][0]
    R2 = data[1][0]
else:
    R1= " ".join(SAMPLES[0])
    R2= " ".join(SAMPLES[1])



shell.prefix("set -eo pipefail; ")

rule all:
    input: 
        f"results/{ORGANISM}_hi-c_library_complexity.txt", 
        f"results/{ORGANISM}_stats_library.txt", 
        f"workflow/samtools/{ORGANISM}.sorted.bam"


rule library_complexity:
    input: f"workflow/samtools/{ORGANISM}.sorted.bam"
    output: f"results/{ORGANISM}_hi-c_library_complexity.txt"
    log:
        os.path.join(output, "logs/library_complexity.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell: 
        "(set +u"
        + " && source {params.source}"
        + " && preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8" 
        + " -seg_len 1000000000 -output {output} {input} ) > {log} 2>&1"

rule get_stats:
    input: "workflow/pairtools/stats.txt"
    output: f"results/{ORGANISM}_stats_library.txt"
    log:
        os.path.join(output, "logs/get_stats.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell: 
        "(set +u" 
        + " && source {params.source}"
        + " && get_qc.py -p {input} > {output} ) > {log} 2>&1"


rule sort_bam:
    input: f"workflow/pairtools/{ORGANISM}.bam"
    output: 
        sort = f"workflow/samtools/{ORGANISM}.sorted.bam", 
        index = f"workflow/samtools/{ORGANISM}.sorted.bam.bai"
    log:
        os.path.join(output, "logs/sort_bam.log")
    threads:
        int(HPC_CONFIG.get_cores("sort_bam"))
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u"
        + " && source {params.source}"
        + " && samtools sort -@{threads} -T /tmp/tmp.bam -o {output.sort} {input}"
        + " && samtools index {output.sort} ) > {log} 2&1"

rule unsorted_bam:
    input: "workflow/pairtools/dedup.pairsam"
    output: 
        map = "workflow/pairtools/mapped.pairs", 
        bam = f"workflow/pairtools/{ORGANISM}.bam"
    log:
        os.path.join(output, "logs/unsorted_bam.log")
    threads:
        int(HPC_CONFIG.get_cores("unsorted_bam")) 
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u"
        + " && source {params.source}" 
        + " && pairtools split --nproc-in {threads} --nproc-out {threads}"
        + " --output-pairs {output.map} --output-sam {output.bam} {input} ) {log} > 2&1"

rule pairtools_dedup:
    input: "workflow/pairtools/sorted.pairsam"
    output: 
        out = "workflow/pairtools/dedup.pairsam", 
        stats =  "workflow/pairtools/stats.txt"
    log:
        os.path.join(output, "logs/pairtools_dedup.log")
    params:
        source = config["source"]["omni-c"]
    threads:
        int(HPC_CONFIG.get_cores("pairtools_dedup"))
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " &&  pairtools dedup --nproc-in {threads} --nproc-out {threads}"
        + " -mark-dups --output-stats {output.stats} --output {output.out} {input}"
        + ") > {log} 2>&1"
    

rule pairtools_sort:
    input:"workflow/pairtools/parsed.pairsam"
    output:"workflow/pairtools/sorted.pairsam"
    params:
        source = config["source"]["omni-c"]
    threads:
        int(HPC_CONFIG.get_cores("pairtools_sort"))
    log:
        os.path.join(output, "logs/pairtools_sort.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " &&  pairtools sort --nproc {threads} --tmpdir=/tmp {input} > {output}"
        + ") > {log} 2>&1"
    

rule unique_unique:
    input: 
        align = f"workflow/bwa/{ORGANISM}_mapped_reads.sam",
        reference = f"{REFERENCE}"
    output: "workflow/pairtools/parsed.pairsam"
    threads:
        HPC_CONFIG.get_cores("unique_unique")
    log:
        os.path.join(output, "logs/unique_unique.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    params:
        source = config["source"]["omni-c"]
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " &&  pairtools parse --min-mapq 40 --walks-policy 5unique" 
        + " --max-inter-align-gap 30"
        + " --nproc-in {threads} --chroms-path {input.reference} {inputs.align}"
        + ") > {log} 2>&1"


rule mapping_to_reference:
    input:
        R1=R1,
        R2=R2,
        reference= f"{REFERENCE}",
        index= f"reference/genome/{ORGANISM}.fasta.fai"
    output: f"workflow/bwa/{ORGANISM}_mapped_reads.sam"
    params:
        source = config["source"]["omni-c"]
    threads:
        int(HPC_CONFIG.get_cores("mapping_to_reference"))
    log:
        os.path.join(output, "logs/mapping_to_reference.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " && bwa mem -5SP -T0 -t {threads} {input.reference} <(zcat {input.R1})"
        + " <(zcat {input.R2}) -o {output}"
        + ") > {log} 2>&1"

rule build_index:
    input: 
        reference = f"{REFERENCE}" , 
        index = f"reference/genome/{ORGANISM}.fasta.fai", 
        genome = f"reference/genome/{ORGANISM}.genome"
    output: f"reference/genome/{ORGANISM}.fasta.amb"
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    log:
        os.path.join(output, "logs/build_index.log")
    shell: 
        "(set +u" 
        + " && source {params.source}" 
        + " && bwa index {input.reference}"
        + ") > {log} 2>&1"

rule build_genome:
    input: f"reference/genome/{ORGANISM}.fasta.fai"
    output: f"reference/genome/{ORGANISM}.genome",
    log:
        os.path.join(output, "logs/build_genome.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u"  
        + " && cut -f1,2 {input} > {output}"
        + ") > {log} 2>&1"
    
rule index_reference:
    input: f"{REFERENCE}"      
    output: f"reference/genome/{ORGANISM}.fasta.fai"   
    log:
        os.path.join(output, "logs/index_reference.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    params:
        source = config["source"]["omni-c"]
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " && samtools faidx {input}"
        + ") > {log} 2>&1"