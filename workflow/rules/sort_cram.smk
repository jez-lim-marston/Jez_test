#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule samtools_sort_cram:
    input:
        bam = "results/align_multi/{sample_id}/Aligned.out.bam",
        ref = config['genome_fasta']        
    output:
        "results/align_multi/{sample_id}/Aligned.sortedByCoord.out.bam",
        "results/align_multi/{sample_id}/Aligned.sortedByCoord.out.cram"
    conda:
        "../envs/utils.yaml"
    threads: 8
    shell:
        '''
tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sample_id}.XXXXXX)
samtools sort -@ {threads} -T $tdir {input.bam} -o {output[0]}
samtools view -C -T {input.ref} -o {output[1]} {output[0]}
samtools index {output[1]}
        '''
