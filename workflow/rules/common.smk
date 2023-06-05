#! /usr/bin/env python
# -*- coding: utf-8 -*-

#read in the metadata, and pharse into sample 
gdc_file = pd.read_csv(config['tcga_metadata'],header = 0)
gdc_file.set_index(config['tcga_sample_id'],inplace = True, drop = False)
tcga_samples = gdc_file[config['tcga_sample_id']].unique().tolist()

""" All samples """
samples = pd.DataFrame(
        {'sample_id': tcga_samples}).set_index("sample_id", drop=False)



