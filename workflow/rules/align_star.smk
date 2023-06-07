#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule star_align_multi:
    input:
        ubam = 'results/ubam/{sample_id}.bam',
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

picard SamToFastq\
  -I {input.ubam}\
  -F $tdir/R1.fq\
  -F2 $tdir/R2.fq\
  --CLIPPING_ATTRIBUTE XT\
  --CLIPPING_ACTION N\
  -NON_PF true     

STAR\
  --runThreadN {threads}\
  --genomeDir {input.genomeDir}\
  --readFilesIn $tdir/R1.fq $tdir/R2.fq\
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

rule samtools_bam_cram:
    input:
        bam = "results/align_multi/{sample_id}/Aligned.out.bam",
        ref = config['genome_fasta']
    output:
        "results/align_multi/{sample_id}/Aligned.out.cram",
    conda:
        "../envs/utils.yaml"
    threads: 2
    shell:
        '''
samtools view -C -T {input.ref} {input.bam} > {output[0]}
        '''
#rule star_align_multi:
#    input:
#        read_1 = 'runs/{sample_id}_1.fastq.gz',
#        read_2 = 'runs/{sample_id}_2.fastq.gz',
#        genomeDir = config['star_align_multi']['genomeDir']
#    output:
#        "results/align_multi/{sample_id}/Aligned.out.bam",
#        "results/align_multi/{sample_id}/ReadsPerGene.out.tab",
#        "results/align_multi/{sample_id}/Log.final.out",
#        "results/align_multi/{sample_id}/SJ.out.tab"
#    conda:
#        "../envs/star.yaml"
#    threads: snakemake.utils.available_cpu_count()
#    shell:
#        '''
#tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sample_id}.XXXXXX)

#STAR\
#  --runThreadN {threads}\
#  --genomeDir {input.genomeDir}\
#  --readFilesIn {input.read_1} {input.read_2}\
#  --readFilesCommand zcat\
#  --outSAMattributes NH HI NM MD AS XS\
#  --outSAMtype BAM Unsorted\
#  --outFileNamePrefix $(dirname {output[0]})/\
#  --quantMode GeneCounts\
#  --outSAMstrandField intronMotif\
#  --outFilterMultimapNmax {config[star_align_multi][outFilterMultimapNmax]}\
#  --winAnchorMultimapNmax {config[star_align_multi][winAnchorMultimapNmax]}\
#  --outSAMunmapped Within KeepPairs
#rm -rf $tdir
#        '''



