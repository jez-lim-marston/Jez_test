#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import csv
import pandas as pd

rule get_sams:
    conda: "../envs/dump.yaml"
    params:
        cmd = lambda wc: allen_readtable.loc[wc.sample_id][config['allen_sample_id']],
        ngc = config['ngc_file']
    output: 
        org = "runs/{umid}.org.bam"
    wildcard_constraints:
        umid = '(UMM|BSSR).+'
    resources:
        mem_mb = config['fastq_mem_mb']
    shell:
        '''
mkdir -p runs/{wildcards.umid}
prefetch {params.cmd} --ngc {params.ngc} -O runs/{wildcards.umid} -X 9999999999999
sam-dump --unaligned runs/{wildcards.umid}/{params.cmd}/{params.cmd}.sra --ngc {params.ngc} | samtools view -bS > {output.org}
        '''

# Default SAM attributes cleared by RevertSam
attr_revertsam = ['NM', 'UQ', 'PG', 'MD', 'MQ', 'SA', 'MC', 'AS']
# SAM attributes output by STAR
attr_star = ['NH', 'HI', 'NM', 'MD', 'AS', 'nM', 'jM', 'jI', 'XS', 'uT']
# Additional attributes to clear
ALN_ATTRIBUTES = list(set(attr_star) - set(attr_revertsam))

rule revert_and_mark_adapters:
    """ Create unmapped BAM (uBAM) from aligned BAM
    """
    input:
        "results/original_bam/{sample_id}.bam"
    output:
        "results/ubam/{sample_id}.bam"
    wildcard_constraints:
        sample_id =  "(SRR)[0-9]+"
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

localrules: make_allen_ubams
rule make_allen_ubams:
    input:
        expand("results/ubam/{sample_id}.bam", sample_id=allen_SAMPLES)
