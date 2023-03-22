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
LONG_READS = list_basenames(config["input_samples"]["long_reads"])
logs = config["logs"]

threading = ''
if config['bwa-mem2'] == 'True':
    bwa = 'bwa-mem2'
    threading = "-p {threads}"
elif config['bwa-mem2'] == 'False':
    bwa = 'bwa'
    
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
        f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam",
        f"{OUTPUT}/workflow/bwa/{ORGANISM}_mapped_reads.sort.bam",
        f"{OUTPUT}/workflow/pretext/{ORGANISM}_unique_mapping.pretext",
        f"{OUTPUT}/workflow/pretext/{ORGANISM}_multi_mapping.pretext",
        f"{OUTPUT}/workflow/cooler/unique_1kb.mcool",
        f"{OUTPUT}/workflow/cooler/all_1kb.mcool",
        f"{OUTPUT}/workflow/tracks/gaps_{ORGANISM}.bedgraph",
        f"{OUTPUT}/workflow/tracks/telomeres_{ORGANISM}.bedgraph",
        f"{OUTPUT}/workflow/tracks/coverage_{ORGANISM}.bedgraph"


rule mosdepth_coverage:
    input:
        f"{OUTPUT}/workflow/samtools/{ORGANISM}_long_reads.sort.bam"
    output:
        f"{OUTPUT}/workflow/tracks/{ORGANISM}.mosdepth.summary.txt",
        bedgz = temp(f"{OUTPUT}/workflow/tracks/{ORGANISM}.per-base.bed.gz"),
        bedgraph = f"{OUTPUT}/workflow/tracks/coverage_{ORGANISM}.bedgraph"
    params:
        source = config["source"]["mosdepth"]
    log:
        os.path.join(logs, "mosdepth_coverage.log")
    threads:
        int(HPC_CONFIG.get_cores("mosdepth_coverage"))
    resources:
        mem_mb = HPC_CONFIG.get_memory("mosdepth_coverage")
    shell:
        "(set + u"
        + " && source {params.source}"
        + f" && cd {OUTPUT}/workflow/tracks"
        + f" && mosdepth -t {{threads}} -b 10000 -x  -m {ORGANISM} {{input}}"
        + " && zcat {output.bedgz} > {output.bedgraph}"
        + " ) > {log} 2>&1"

rule index_bam:
    input:
        f"{OUTPUT}/workflow/samtools/{ORGANISM}_long_reads.sort.bam"
    output:
        f"{OUTPUT}/workflow/samtools/{ORGANISM}_long_reads.sort.bam.bai"
    params:
        source = config["source"]["omni-c"]
    log:
        os.path.join(logs, "index_bam.log")
    threads:
        int(HPC_CONFIG.get_cores("index_bam"))
    resources:
        mem_mb = HPC_CONFIG.get_memory("index_bam")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && samtools index -@ {threads} {input}"
        + " ) > {log} 2>&1"

rule merge_bam:
    input:
        expand(f"{OUTPUT}/workflow/samtools/{{long_reads}}.sort.bam", long_reads = LONG_READS)
    output:
        bam = f"{OUTPUT}/workflow/samtools/{ORGANISM}_long_reads.sort.bam",
        index = f"{OUTPUT}/workflow/samtools/{ORGANISM}_long_reads.sort.bam"
    log:
        os.path.join(logs, "merge_bam.log")
    params:
        source = config["source"]["omni-c"]
    threads:
        int(HPC_CONFIG.get_cores("merge_bam"))
    resources:
        mem_mb = HPC_CONFIG.get_memory("merge_bam")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && samtools merge -@ {threads} {output.bam} {input}"
        + f" && samtools index -@ {{threads}} {OUTPUT}/workflow/samtools/{ORGANISM}_long_reads.sort.bam"
        + " ) > {log} 2>&1"

rule minimap2:
    input:
        reads = expand(f"{OUTPUT}/reads/{{long_reads}}", long_reads=LONG_READS),
        reference = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta",
    output:
        expand(f"{OUTPUT}/workflow/samtools/{{long_reads}}.sort.bam", long_reads=LONG_READS),
    log:
        os.path.join(logs, "minimap2.log")
    params:
        source = config["source"]["minimap2"],
        source_2 = config["source"]["omni-c"],
    threads:
        int(HPC_CONFIG.get_cores("minimap2"))
    resources:
        mem_mb = HPC_CONFIG.get_memory("minimap2")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && source {params.source_2}"
        + " && minimap2 -t {threads} -ax map-hifi {input.reference} {input.reads} " 
        + " | samtools view -@ {threads} -bS - | samtools sort -@ {threads} -m 4G > "
        + " {output}"
        + " ) > {log} 2>&1"

rule telomeres_reformatting:
    input:
        f"{OUTPUT}/workflow/telomeres/{ORGANISM}_telomeric_repeat_windows.tsv"
    output:
        f"{OUTPUT}/workflow/tracks/telomeres_{ORGANISM}.bedgraph"
    shell:""" 
        awk '$6=$2-10000  {{print $1,$6,$2,$3}}' {input} | sed "s/\ /\\t/g" > {output}
        """

rule telomeres:
    input:
        fasta = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta",
    output:
        f"{OUTPUT}/workflow/telomeres/{ORGANISM}_telomeric_repeat_windows.tsv"
    log:
        os.path.join(logs, "telomeres.log")
    params:
        source = config["source"]["tidk"],
        prefix = f"{ORGANISM}",
        #sequence = config["input_samples"]["telomeres"],
        dir = f"{OUTPUT}/workflow/telomeres"
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && tidk search --extension tsv -s TTAGGG -o {params.prefix} "
        + " --dir {params.dir}  {input.fasta}"
        + " && sed  -i 1d {output} " 
        + " ) > {log} 2>&1"


rule gaps_reformatting:
    input:
        f"{OUTPUT}/workflow/tracks/gaps.bed3"
    output:
        f"{OUTPUT}/workflow/tracks/gaps_{ORGANISM}.bedgraph"
    shell: """ 
        awk ' $4=$3-$2 {{ print $0}}' {input} > {output} 
        """

rule gaps:
    input:
        fasta = f"{OUTPUT}/reference/genome/{ORGANISM}.fasta",
    output:
        temp(f"{OUTPUT}/workflow/tracks/gaps.bed3")
    log:
        os.path.join(logs, "gaps.log")
    params:
        source = config["source"]["seqtk"]
    shell:
        "(set +u" 
        + " && source {params.source}"
        + f" && mkdir -p {OUTPUT}/workflow/tracks"
        + " && seqtk cutN -g -n 10 {input.fasta} > {output} "
        + " ) > {log} 2>&1"


rule uniquemapping_pretext:
    input:
        f"{OUTPUT}/workflow/samtools/{ORGANISM}.sorted.bam"
    output: 
        f"{OUTPUT}/workflow/pretext/{ORGANISM}_unique_mapping.pretext"
    params: 
        source = config["source"]["omni-c"],
        source_2 = config["source"]["pretext"]
    threads:
        int(HPC_CONFIG.get_cores("uniquemapping_pretext"))
    log:
        os.path.join(logs, "uniquemapping_pretext.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + " && source {params.source_2}"
        + " && samtools view -h {input} | PretextMap --mapq 0 -o {output}"
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
        mem_mb = HPC_CONFIG.get_memory("sort_bam")
    shell:
        "(set +u"
        + " && source {params.source}"
        + f" && samtools sort -m 4G -@ {{threads}} -T {OUTPUT}/workflow/tmp.bam"
        + " -o {output.sort} {input}"
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
        mem_mb = HPC_CONFIG.get_memory("unsorted_bam")
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
    input:
        f"{OUTPUT}/workflow/pairtools/unique.pairs.gz"
    output:
        f"{OUTPUT}/workflow/pairtools/sorted.pairbam"
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
    
rule cooler_unique:
    input:
        pairs = f"{OUTPUT}/workflow/pairtools/unique.sort.pairs.gz",
        index = f"{OUTPUT}/workflow/pairtools/unique.sort.pairs.gz.px2",
        genome = f"{OUTPUT}/reference/genome/{ORGANISM}.genome"
    output:
        cool = temp(f"{OUTPUT}/workflow/cooler/unique_1kb.cool"),
        mcool = f"{OUTPUT}/workflow/cooler/unique_1kb.mcool"
    threads:
        HPC_CONFIG.get_cores("cooler_unique")
    log:
        os.path.join(logs, "cooler_unique.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("cooler_unique")
    params: 
        source = config["source"]["cooler"],
    shell:
        "(set +u" 
        + " && source {params.source}" 
        + f" && cooler cload pairix --assembly {ORGANISM}"
        + " -p {threads} {input.genome}:1000 {input.pairs} {output.cool}"
        + f" && cooler zoomify --balance -r 1000N -p {{threads}} {OUTPUT}/workflow/cooler/unique_1kb.cool"
        + " ) > {log} 2>&1"
        

rule cooler_all:
    input:
        pairs = f"{OUTPUT}/workflow/pairtools/all.sort.pairs.gz",
        genome = f"{OUTPUT}/reference/genome/{ORGANISM}.genome",
        index = f"{OUTPUT}/workflow/pairtools/all.sort.pairs.gz.px2",
    output:
        cool = temp(f"{OUTPUT}/workflow/cooler/all_1kb.cool"),
        mcool = f"{OUTPUT}/workflow/cooler/all_1kb.mcool"
    threads:
        HPC_CONFIG.get_cores("cooler_all")
    log:
        os.path.join(logs, "unique_ucooler_allnique.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("cooler_all")
    params: 
        source = config["source"]["cooler"],
    shell:
        "(set +u" 
        + " && source {params.source}"
        + f" && cooler cload pairix --assembly {ORGANISM}"
        + " -p {threads} {input.genome}:1000 {input.pairs} {output.cool}" 
        + f" && cooler zoomify --balance -r 1000N -p {{threads}} {OUTPUT}/workflow/cooler/all_1kb.cool"
        + " ) > {log} 2>&1"
        


rule unique_unique:
    input: 
        align = f"{OUTPUT}/workflow/bwa/{ORGANISM}_mapped_reads.sort.bam",
        reference = f"{OUTPUT}/reference/genome/{ORGANISM}.genome"
    output: 
        unique = temp(f"{OUTPUT}/workflow/pairtools/unique.pairs"),
        compress = f"{OUTPUT}/workflow/pairtools/unique.sort.pairs.gz",
        stats = f"{OUTPUT}/workflow/pairtools/unique.pairs.stats",
        index = f"{OUTPUT}/workflow/pairtools/unique.sort.pairs.gz.px2"
    threads:
        HPC_CONFIG.get_cores("unique_unique")
    log:
        os.path.join(logs, "unique_unique.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("mapping_to_reference")
    params: 
        source = config["source"]["omni-c"],
        source2 = config["source"]["cooler"]
    shell:
        "(set +u"
        + f" && mkdir -p {OUTPUT}/workflow/pairtools/" 
        + " && source {params.source}" 
        + " &&  pairtools parse --min-mapq 40 --walks-policy 5unique" 
        + " --max-inter-align-gap 30 -o {output.unique}"
        + " --output-stats {output.stats}"
        + " --nproc-in {threads} --chroms-path {input.reference} {input.align}"
        + f" &&  pairtools sort --nproc {{threads}} --tmpdir={OUTPUT}/tmp "
        + f" {OUTPUT}/workflow/pairtools/unique.pairs > {OUTPUT}/workflow/pairtools/unique.sort.pairs"
        + " && source {params.source2}"
        + f" && bgzip {OUTPUT}/workflow/pairtools/unique.sort.pairs"
        + f" && pairix -p pairs -f {OUTPUT}/workflow/pairtools/unique.sort.pairs.gz"
        + " || true"
        + " ) > {log} 2>&1"



rule multimapping_pairtools:
    input: 
        align = f"{OUTPUT}/workflow/bwa/{ORGANISM}_mapped_reads.sort.bam",
        reference = f"{OUTPUT}/reference/genome/{ORGANISM}.genome"
    output: 
        unique = temp(f"{OUTPUT}/workflow/pairtools/all.pairs"),
        compress = f"{OUTPUT}/workflow/pairtools/all.sort.pairs.gz",
        stats = f"{OUTPUT}/workflow/pairtools/all.pairs.stats",
        index = f"{OUTPUT}/workflow/pairtools/all.sort.pairs.gz.px2"
    threads:
        HPC_CONFIG.get_cores("multimapping_pairtools")
    log:
        os.path.join(logs, "multimapping_pairtools.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("multimapping_pairtools")
    params: 
        source = config["source"]["omni-c"],
        source2 = config["source"]["cooler"]
    shell:
        "(set +u"
        + f" && mkdir -p {OUTPUT}/workflow/pairtools/" 
        + " && source {params.source}" 
        + " &&  pairtools parse --min-mapq 0 --walks-policy all" 
        + " --max-inter-align-gap 30 -o {output.unique}"
        + " --output-stats {output.stats}"
        + " --nproc-in {threads} --chroms-path {input.reference} {input.align}"
        + f" &&  pairtools sort --nproc {{threads}} --tmpdir={OUTPUT}/tmp"
        + f" {OUTPUT}/workflow/pairtools/all.pairs > {OUTPUT}/workflow/pairtools/all.sort.pairs"
        + " && source {params.source2
        + f" && bgzip {OUTPUT}/workflow/pairtools/all.sort.pairs"
        + f" && pairix -p pairs -f {OUTPUT}/workflow/pairtools/all.sort.pairs.gz"
        + " || true "
        + " ) > {log} 2>&1"


rule multimapping_pretext:
    input:
        f"{OUTPUT}/workflow/bwa/{ORGANISM}_mapped_reads.sort.bam"
    output: 
        f"{OUTPUT}/workflow/pretext/{ORGANISM}_multi_mapping.pretext"
    params: 
        source = config["source"]["omni-c"],
        source_2 = config["source"]["pretext"],
    threads:
        int(HPC_CONFIG.get_cores("multimapping_pretext"))
    log:
        os.path.join(logs, "multimapping_pretext.log")
    resources:
        mem_mb = HPC_CONFIG.get_memory("multimapping_pretext")
    shell:
        "(set +u" 
        + " && source {params.source}"
        + " && source {params.source_2}"
        + " && samtools view -h {input} | PretextMap -o {output}"
        + " ) > {log} 2>&1"


rule mapping_to_reference:
    input:
        R1 = expand(f"{OUTPUT}/reads/{r1}", r1=R1),
        R2 = expand(f"{OUTPUT}/reads/{r2}", r2=R2),
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
        + " && {params.bwa} index {params.threading}"
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

