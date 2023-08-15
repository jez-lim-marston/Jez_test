#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule telescope:
    input:
        mbam = "results/align_multi/{sample_id}/Aligned.out.bam",
        annot = config['telescope']['annotation']
    output:
        "results/telescope/{sample_id}/report.tsv",
        "results/telescope/{sample_id}/updated.bam",
        "results/telescope/{sample_id}/other.bam"
    log:
        "results/telescope/{sample_id}/telescope.log"        
    conda:
        "../envs/telescope.v1.yaml"
    shell:
        '''
tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sample_id}.XXXXXX)

telescope assign\
 --exp_tag inform\
 --theta_prior 200000\
 --max_iter 1000\
 --updated_sam\
 --outdir $tdir\
 {input.mbam}\
 {input.annot}\
 2>&1 | tee {log[0]}

mkdir -p $(dirname {output[0]})
mv $tdir/inform-telescope_report.tsv {output[0]}
mv $tdir/inform-updated.bam {output[1]}
mv $tdir/inform-other.bam {output[2]}

chmod 600 {output[1]}
rm -rf $tdir
        '''
