#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import csv
import pandas as pd

#download the unaligned bams
#rule get_bams:
#    conda: "../envs/dump.yaml"
#    params:
#        SRR = lambda wc: allen_readtable.loc[wc.sample_id][config['allen_sample_id']],
#        ngc = config['ngc_file'],
#    output: 
#        org = temp("runs/{sample_id}.org.bam")
#    wildcard_constraints:
#        sample_id = "(SRR)[0-9]+"
#    resources:
#        mem_mb = config['fastq_mem_mb']
#    shell:
#        '''
#mkdir -p runs/{wildcards.sample_id}
#prefetch {params.SRR} --ngc {params.ngc} -O runs/{wildcards.sample_id} -X 9999999999999
#sam-dump --unaligned runs/{wildcards.sample_id}/{params.SRR}/{params.SRR}.sra --ngc {params.ngc} | samtools view -bS > {output.org}
#        '''

# get fastq for hugo samples, all reads are paired
rule get_bams:
    conda: "../envs/dump.yaml"
    params:
        SRR = lambda wc: hugo_readtable.loc[wc.sample_id][config['hugo_sample_id']],
        ngc = config['ngc_file'],
    output: 
        read_1 = "runs/{sample_id}_1.fastq.gz",
        read_2 = "runs/{sample_id}_2.fastq.gz",
        sra = temp("runs/{wildcards.sample_id}/{params.SRR}/{params.SRR}.sra ")
    wildcard_constraints:
        sample_id = "(SRR)[0-9]+"
    resources:
        mem_mb = config['fastq_mem_mb']
    shell:
        '''
mkdir -p runs/{wildcards.sample_id}
prefetch {params.SRR} --ngc {params.ngc} -O runs/{wildcards.sample_id} -X 9999999999999
fastq-dump --gzip --split-3 -O runs --ngc {params.ngc} runs/{wildcards.sample_id}/{params.SRR}/{params.SRR}.sra 
        '''

#localrules: allen_fastq
#rule allen_fastq:
#    input:
#        expand("runs/{sample_id}_2.fastq.gz", sample_id=allen_SAMPLES)

localrules: hugo_fastq
rule hugo_fastq:
    input:
        expand("runs/{sample_id}_2.fastq.gz", sample_id=hugo_SAMPLES)