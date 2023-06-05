#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import csv
import pandas as pd

def pref_command(wc):
    SRR = METADATA.loc[wc.umid][config['SOURCE']]
    return SRR

#The files were submitted as bam files.
#The are unorganized - the header should be containing CB/UB information
rule get_sams:
    conda: "../envs/fastq.yaml"
    params:
        cmd = pref_command,
        ngc = config['ngc_file'],
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
        sample_id = "TCGA\\-..\\-[A-Z]...\\-..[A-Z]"
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

localrules: make_tcga_ubams
rule make_tcga_ubams:
    input:
        expand("results/ubam/{sample_id}.bam", sample_id=tcga_samples)