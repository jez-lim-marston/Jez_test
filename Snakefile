#! /usr/bin/env python
# -*- coding: utf-8 -*-
from snakemake.utils import min_version
import os
import pandas as pd
import re

min_version("6.4.1")

configfile: "config/config.yaml"

""" Allen samples """
allen_readtable = pd.read_csv(config['allen_metadata'], header = 0)
allen_readtable.set_index(config['allen_SRR'], inplace = True, drop = False)
allen_SAMPLES = allen_readtable[config['allen_sample_id']].unique().tolist()
allen_readtable.set_index(config['allen_sample_id'], inplace = True, drop = False)

""" Hugo samples """
hugo_readtable = pd.read_csv(config['hugo_metadata'], header = 0)
hugo_readtable.set_index(config['hugo_SRR'], inplace = True, drop = False)
hugo_SAMPLES = allen_readtable[config['hugo_sample_id']].unique().tolist()
hugo_readtable.set_index(config['hugo_sample_id'], inplace = True, drop = False)

include: "workflow/rules/download_and_make_ubam.smk"
include: "workflow/rules/align_star.smk"
include: "workflow/rules/sort_cram.smk"
include: "workflow/rules/telescope.smk"
include: "workflow/rules/stringtie.smk"

localrules: allen_sample_download
rule allen_sample_download:
    input:
        expand("results/ubam/{sample_id}.bam", sample_id=allen_SAMPLES)

localrules: hugo_sample_download
rule hugo_sample_download:
    input:
        expand("results/ubam/{sample_id}.bam", sample_id=hugo_SAMPLES)

## when job completed, remove the ubam to save space
localrules: star_align
rule star_align:
    input:
        "results/align_multi/{sample_id}/Aligned.out.bam",
        "results/align_multi/{sample_id}/ReadsPerGene.out.tab",
        "results/align_multi/{sample_id}/Aligned.sortedByCoord.out.cram",
        "results/ubam/{sample_id}.bam",
    output:
        touch("results/complete/{sample_id}_star_align.txt")
    shell:
        '''
rm {input[3]}
        '''               

localrules: star_complete_allen
rule star_complete_allen:
    input:
        expand("results/complete/{s}_star_align.txt", s=allen_SAMPLES)

localrules: star_complete_hugo
rule star_complete_hugo:
    input:
        expand("results/complete/{s}_star_align.txt", s=hugo_SAMPLES)

localrules: sample_complete
rule sample_complete:
    input:
        "results/align_multi/{sample_id}/Aligned.out.bam",
        rules.telescope.output,
        rules.stringtie.output
    output:
        touch("results/complete/{sample_id}.txt")
    shell:
        '''
rm {input[0]}
        '''      

localrules: allen_all_complete
rule allen_all_complete:
    input:
        expand("results/complete/{s}.txt", s=allen_SAMPLES)

localrules: hugo_all_complete
rule hugo_all_complete:
    input:
        expand("results/complete/{s}.txt", s=hugo_SAMPLES)
