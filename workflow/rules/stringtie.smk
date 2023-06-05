#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule stringtie:
    input:
        bam = "results/align_multi/{sample_id}/Aligned.sortedByCoord.out.bam",
        ref = config['genome_fasta'],
        annot = config['annotation_gtf']
    output:
        "results/stringtie/{sample_id}/transcripts.gtf"
    conda: 
        "../envs/stringtie.yaml"
    threads: 2
    shell:
        '''
stringtie\
 -p {threads}\
 -u \
 -c 2.5\
 -s 2.5\
 -j 2\
 -f 0.05\
 -M 1\
 -G {input.annot}\
 -o {output[0]}\
 {input.bam}
        '''
