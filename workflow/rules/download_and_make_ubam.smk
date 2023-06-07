#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import csv
import pandas as pd

#download data from van allen dataset, all samples paired
rule get_bams:
    conda: "../envs/dump.yaml"
    params:
        SRR = lambda wc: allen_readtable.loc[wc.sample_id][config['allen_sample_id']],
        ngc = config['ngc_file'],
    output: 
        org = temp("runs/{sample_id}.org.bam")
    wildcard_constraints:
        sample_id = "(SRR)[0-9]+"
    resources:
        mem_mb = config['fastq_mem_mb']
    shell:
        '''
mkdir -p runs/{wildcards.sample_id}
prefetch {params.SRR} --ngc {params.ngc} -O runs/{wildcards.sample_id} -X 9999999999999
sam-dump --unaligned runs/{wildcards.sample_id}/{params.SRR}/{params.SRR}.sra --ngc {params.ngc} | samtools view -bS > {output.org}
        '''

# Default SAM attributes cleared by RevertSam
attr_revertsam = ['NM', 'UQ', 'PG', 'MD', 'MQ', 'SA', 'MC', 'AS']
# SAM attributes output by STAR
attr_star = ['NH', 'HI', 'NM', 'MD', 'AS', 'nM', 'jM', 'jI', 'XS', 'uT']
# Additional attributes to clear
ALN_ATTRIBUTES = list(set(attr_star) - set(attr_revertsam))

##Create unmapped BAM (uBAM) from aligned BAM
rule revert_and_mark_adapters:
    input:
        "runs/{sample_id}.org.bam"
    output:
        "results/ubam/{sample_id}.bam"
    wildcard_constraints:
        sample_id = "(SRR)[0-9]+"
    conda:
        "../envs/utils.yaml"
    params:
        attr_to_clear = expand("--ATTRIBUTE_TO_CLEAR {a}", a=ALN_ATTRIBUTES),
        tmpdir = config['tmp']
    shell:
        '''
picard RevertSam\
 -I {input[0]}\
 -O {output[0]}\
 --SANITIZE true\
 --COMPRESSION_LEVEL 0\
 {params.attr_to_clear}\
 --TMP_DIR {params.tmpdir}
chmod 600 {output[0]}
        '''

localrules: allen_bam
rule allen_bam:
    input:
        expand("results/ubam/{sample_id}.bam", sample_id=allen_SAMPLES)

# get fastq for hugo samples, all reads are paired
#rule get_bams:
#    conda: "../envs/dump.yaml"
#    params:
#        SRR = lambda wc: hugo_readtable.loc[wc.sample_id][config['hugo_sample_id']],
#        ngc = config['ngc_file'],
#    output:
#        read_1 = "runs/{sample_id}_1.fastq.gz",
#        read_2 = "runs/{sample_id}_2.fastq.gz"
#    wildcard_constraints:
#        sample_id = "(SRR)[0-9]+"
#    resources:
#        mem_mb = config['fastq_mem_mb']
#    shell:
#        '''
#mkdir -p runs/{wildcards.sample_id}
#prefetch {params.SRR} --ngc {params.ngc} -O runs/{wildcards.sample_id} -X 9999999999999
#fastq-dump --gzip --split-3 -O runs --ngc {params.ngc} runs/{wildcards.sample_id}/{params.SRR}/{params.SRR}.sra
#        '''

#localrules: hugo_fastq
#rule hugo_fastq:
#    input:
#        expand("runs/{sample_id}_2.fastq.gz", sample_id=hugo_SAMPLES)