import os
import sys
import json
import csv
from snakemake.utils import min_version
min_version("7.10")

from eihic.scripts.hpc_config import HpcConfig
HPC_CONFIG = HpcConfig(config["hpc_config"])

# declare variables
jira_id = config["jira"]["jira_id"] # command line
NOTIFY = not config["notify"] # command line
if not Path(config["jira"]["password_file"]).is_file() or not jira_id:
    NOTIFY = False
output = os.path.abspath(config["output"])
logs_dir = os.path.abspath(config["logs"])
output = os.path.abspath(config["output"]
#save the data to a list: R1 forward reads, R2 reverse reads, and sample name/ organism name
R1= config["R1"]
R2= config["R2"]
REFERENCE= config["reference"]
ORGANISM= config["organism"]


def get_read_list(data):
#create a list of reads for the config file of HiCUP
    global reads_list
    reads_list = ""
    for i in range(len(R1)):

        reads_list = reads_list + output + '/references/reads/' + R1[i] + "\n" + output + '/references/reads/' + R2[i] + "\n\n"
   
    return reads_list

shell.prefix("set -eo pipefail; ")


rule all:
    input: "done"


rule hicup_run:
    input: 
        config = f"{output}/workflow/hicup/hicup_config.txt",
        index = f"{output}/workflow/bowtie2/{ORGANISM}.1.bt2"
    output: "done"
    log:
        os.path.join(output, "logs/hicup_run/hicup_run.log")
    threads:
        int(HPC_CONFIG.get_cores("hicup_run")) 
    params:
        source = config["source"]["chic"]
    shell:
        "(set +u"
        + " && source {params.source}"
        + " && hicup --config {input.config}"
        + " && touch done ) > {log} 2>&1"

get_read_list(data) # gets the reads formatted to use in the config file 
rule hicup_config:
    input: 
        digest = f"{output}/workflow/hicup/{ORGANISM}_arima_digest.txt"
    output: f"{output}/workflow/hicup/hicup_config.txt"
    params:
        outdir = f"{output}/workflow/hicup/",
        keep = str(config["hicup"]["keep"]),
        quiet = str(config["hicup"]["quiet"]),
        gzip = str(config["hicup"]["compress"]),
        bowtie2 = "/opt/software/bin/bowtie2",
        index = f"{output}/workflow/bowtie2/{ORGANISM}",
        formatting = str(config["hicup"]["formatting"]),
        shortest = str(config["hicup"]["shortest_sequence"]),
        longest = str (config["hicup"]["longest_ditag"]),
        source = config["source"]["chic"]
    threads:
        int(HPC_CONFIG.get_cores("def"))
    log:
        os.path.join(output, "logs/hicup_config/hicup_config.log")
    shell:
        "(set +u"
        + " && source {params.source}"
        + f" && echo 'Outdir: {{params.outdir}}\nThreads: {{threads}}\nQuiet: {{params.quiet}}\nKeep: {{params.quiet}}\nZip: {{params.gzip}}\nBowtie2: {{params.bowtie2}}\nIndex: {{params.index}}\nDigest: {{input}}\nFormat: {{params.formatting}}\nLongest: {{params.longest}}\nShortest: {{params.shortest}}\n\n{reads_list}\n' > {{output}}"
        + " ) > {log} 2>&1"


rule hicup_digester:
    input: f"{output}/{REFERENCE}"
    output:f"{output}/workflow/hicup/{ORGANISM}_arima_digest.txt"
    log:
        os.path.join(output, "logs/hicup_digester/hicup_digester.log")
    threads:
        int(HPC_CONFIG.get_cores("def"))
    params:
        source = config["source"]["chic"]
    shell:
        "(set +u" 
        + " && source {params.source}"
        + f" && hicup_digester --genome {ORGANISM}"
        + f" --arima --outdir {output}/workflow/hicup/"
        + " {input}"
        + f" && mv {output}/workflow/hicup/* {{output}} ) > {{log}} 2>&1"

rule bwt_2_index:
    input: f"{output}/{REFERENCE}"
    output:f"{output}/workflow/bowtie2/{ORGANISM}.1.bt2"
    benchmark:f"benchmark/bowtie2/{ORGANISM}_index.txt"
    log:
        os.path.join(output, "logs/bwt_2_index/bwt_2_index.log")
    threads:
        int(HPC_CONFIG.get_cores("bwt_2_index"))
    params:
        source = config["source"]["chic"]
    shell:
        "(set +u"
        + " && source {params.source}"
        + f" && cd {output}/workflow/bowtie2/"
        + " && bowtie2-build --threads {threads} {input}"
        + f" {ORGANISM} ) > {{log}} 2>&1" 
