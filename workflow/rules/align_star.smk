#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule star_align_multi:
    input:
        read_1 = 'runs/{sample_id}_1.fastq.gz',
        read_2 = 'runs/{sample_id}_2.fastq.gz',
        genomeDir = config['star_align_multi']['genomeDir']
    output:
        "results/align_multi/{sample_id}/Aligned.out.bam",
        "results/align_multi/{sample_id}/ReadsPerGene.out.tab",
        "results/align_multi/{sample_id}/Log.final.out",
        "results/align_multi/{sample_id}/SJ.out.tab"
    conda:
        "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sample_id}.XXXXXX)

STAR\
  --runThreadN {threads}\
  --genomeDir {input.genomeDir}\
  --readFilesIn {input.read_1} {input.read_2}\
  --readFilesCommand zcat\
  --outSAMattributes NH HI NM MD AS XS\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $(dirname {output[0]})/\
  --quantMode GeneCounts\
  --outSAMstrandField intronMotif\
  --outFilterMultimapNmax {config[star_align_multi][outFilterMultimapNmax]}\
  --winAnchorMultimapNmax {config[star_align_multi][winAnchorMultimapNmax]}\
  --outSAMunmapped Within KeepPairs
rm -rf $tdir
        '''



