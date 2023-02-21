#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Snakemake for arima Hi-C libraries QC
"""

import os
import sys
import json
import csv

from eihic.scripts.hpc_config import HpcConfig
HPC_CONFIG = HpcConfig(config["hpc_config"])

# declare variables
jira_id = config["jira"]["jira_id"] # command line
NOTIFY = not config["notify"] # command line
if not Path(config["jira"]["password_file"]).is_file() or not jira_id:
    NOTIFY = False
R1 = config["input_samples"]["R1"]
R2 = config["input_samples"]["R2"]
REFERENCE = config["input_samples"]["reference"]
ORGANISM = config["input_samples"]["organism"]
OUTPUT = os.path.abspath(config["output"])

cwd = os.getcwd()
logs_dir = os.path.abspath(config["logs"])

def get_read_list(R1, R2):
#create a list of reads for the config file of HiCUP
    reads_list = ""
    if len(R1) == 1:
         reads_list =  OUTPUT + '/reads/' + R1[0] + "\n" + OUTPUT + '/reads/' + R2[0] + "\n\n"
    else:
        for i in range(len(R1)):
            reads_list = reads_list + OUTPUT + '/reads/' + R1[i] + "\n" + OUTPUT + '/reads/' + R2[i] + "\n\n"

    return reads_list

shell.prefix("set -eo pipefail; ")


rule all:
    input: "done"

rule hicup_run:
    input: 
        config = f"{OUTPUT}/workflow/hicup/hicup_config.txt",
        index = f"{OUTPUT}/workflow/bowtie2/{ORGANISM}.1.bt2"
    output: "done"
    log:
        os.path.join(OUTPUT, "logs/hicup_run/hicup_run.log")
    threads:
        int(HPC_CONFIG.get_cores("hicup_run")) 
    params:
        source = config["source"]["arima"]
    shell:
        "(set +u"
        + " && source {params.source}"
        + " && hicup --config {input.config}"
        + " && touch done ) > {log} 2>&1"

reads_list = get_read_list(R1,R2) # gets the reads formatted to use in the config file 

rule hicup_config:
    input:
        digest = f"{OUTPUT}/workflow/hicup/{ORGANISM}_arima_digest.txt"
    output: f"{OUTPUT}/workflow/hicup/hicup_config.txt"
    params:
        outdir = f"{OUTPUT}/workflow/hicup/",
        keep = str(config["hicup"]["keep"]),
        quiet = str(config["hicup"]["quiet"]),
        gzip = str(config["hicup"]["compress"]),
        bowtie2 = "/opt/software/bin/bowtie2",
        index = f"{OUTPUT}/workflow/bowtie2/{ORGANISM}",
        formatting = str(config["hicup"]["formatting"]),
        shortest = str(config["hicup"]["shortest_sequence"]),
        longest = str (config["hicup"]["longest_ditag"]),
        source = config["source"]["arima"],
        reads_list = reads_list
    threads:
        int(HPC_CONFIG.get_cores("def"))
    log:
        os.path.join(OUTPUT, "logs/hicup_config/hicup_config.log")
    shell:
        "(set +u"
        + " && source {params.source}"
        + " && echo 'Outdir: {params.outdir}\nThreads: {threads}\nQuiet: {params.quiet}\nKeep: {params.keep}\nZip: {params.gzip}\nBowtie2:  {params.bowtie2}\nIndex: {params.index}\nDigest: {input}\nFormat:  {params.formatting}\nLongest: {params.longest}\nShortest:  {params.shortest}\n\n{params.reads_list}\n' > {output}"
        + " ) > {log} 2>&1"


rule hicup_digester:
    input: f"{OUTPUT}/reference/genome/{REFERENCE}"
    output: f"{OUTPUT}/workflow/hicup/{ORGANISM}_arima_digest.txt"
    log:
        os.path.join(OUTPUT, "logs/hicup_digester/hicup_digester.log")
    threads:
        int(HPC_CONFIG.get_cores("def"))
    params:
        source = config["source"]["arima"]
    shell:
        "(set +u" 
        + " && source {params.source}"
        + f" && hicup_digester --genome {ORGANISM}"
        + f" --arima --outdir {OUTPUT}/workflow/hicup/"
        + " {input}"
        + f" && mv {OUTPUT}/workflow/hicup/* {{output}} ) > {{log}} 2>&1"

rule bwt_2_index:
    input: f"{OUTPUT}/reference/genome/{REFERENCE}"
    output: f"{OUTPUT}/workflow/bowtie2/{ORGANISM}.1.bt2"
    log:
        os.path.join(OUTPUT, "logs/bwt_2_index/bwt_2_index.log")
    threads:
        int(HPC_CONFIG.get_cores("bwt_2_index"))
    params:
        source = config["source"]["arima"]
    shell:
        "(set +u"
        + " && source {params.source}"
        + f" && cd {OUTPUT}/workflow/bowtie2/"
        + " && bowtie2-build --threads {threads} {input}"
        + f" {ORGANISM} ) > {{log}} 2>&1" 
