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

# Save the data to a list: R1 forward reads, R2 reverse reads, and sample name/ organism name
REFERENCE = config["input_samples"]["reference"]
ORGANISM = config["input_samples"]["organism"]
OUTPUT = config["output"]
logs = config["logs"]


rule all:
    input: "done"

rule coverage:
    input:
        bam = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta",
    output:
        f"{OUTPUT}/workflow/yahs/yahs_{ORGANISM}",
    log:
        os.path.join(logs, "scaffold_yahs.log")
    params:
        source = config["source"]["yahs"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("scaffold_yahs")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && yahs --no-mem-check -r 1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000" 
        + " -o {output} {input.fasta} {input.bam}"
        + " ) > {log} 2>&1"

rule telomeres:
    input:
        dir = f"{OUTPUT}/workflow/telomeres",
        fasta = f"{OUTPUT}/reference/yahs/{ORGANISM}_scaffolds_final.fa",
    output:
        f"{OUTPUT}/workflow/telomeres/{ORGANISM}_telomere.track"
    log:
        os.path.join(logs, "scaffold_yahs.log")
    params:
        source = config["source"]["tidk"],
        prefix = f"{OUTPUT}/workflow/telomeres/{ORGANISM}"
    resources:
        mem_mb = HPC_CONFIG.get_memory("telomeres")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && tidk search --extension tsv -s {params.sequence} -o {output} "
        + " --dir {input.dir} {input.fasta} "
        + " -o {output} {input.fasta} {input.bam}"
        + f" && sed  1d {OUTPUT}/workflow/telomeres/{ORGANISM}_telomeric_repeat_windows.tsv | "
        + """ awk "\$6=\$2-10000  {{print \$1,\$6,\$2,\$3}} " | sed "s/\ /\\t/g" > {output}"""
        + " ) > {log} 2>&1"

rule repeats:
    input:
        bam = f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam",
        fasta = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta",
    output:
        f"{OUTPUT}/workflow/yahs/yahs_{ORGANISM}"
    log:
        os.path.join(logs, "scaffold_yahs.log")
    params:
        source = config["source"]["yahs"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("scaffold_yahs")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && yahs --no-mem-check -r 1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000" 
        + " -o {output} {input.fasta} {input.bam}"
        + " ) > {log} 2>&1"

rule gaps:
    input:
        fasta = f"{OUTPUT}/reference/yahs/{ORGANISM}.fasta",
    output:
        f"{OUTPUT}/workflow/yahs/yahs_{ORGANISM}"
    log:
        os.path.join(logs, "scaffold_yahs.log")
    params:
        source = config["source"]["seqtk"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("scaffold_yahs")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && seqtk cutN -g -n 10 {input.fasta} | "
        + """ awk ' \$4=\$3-\$2 {{ print \$0}} > {output} """
        + " ) > {log} 2>&1"

rule scaffold_yahs:
    input:
        bam = f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam",
        fasta = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta",
    output:
        f"{OUTPUT}/workflow/yahs/yahs_{ORGANISM}"
    log:
        os.path.join(logs, "scaffold_yahs.log")
    params:
        source = config["source"]["yahs"]
    resources:
        mem_mb = HPC_CONFIG.get_memory("scaffold_yahs")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && yahs --no-mem-check -r 1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000" 
        + " -o {output} {input.fasta} {input.bam}"
        + " ) > {log} 2>&1"