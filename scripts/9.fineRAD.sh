########################################
# 9.fineRAD                            #
# Author: Flavio Piras                 #
# Reference: ....                      #
#                                      #
########################################



#!/bin/bash

 populations \ #populations is a module in Stacks package
    -P /path/to/gstacks_output \
    -M /path/to/pop_map.txt \ #file with sample_id \t population_locality format 
    -O /path/to/output_folder \
    -t 16 \
    --min-gt-depth 10 \
    -R 0.6 \
    -p 2 \ #minimum number of populations a locus must be present in to process a locus.
    --max-obs-het 0.7 \
    --vcf-all #for pixy analysis which requires not just the SNPs 

RADpainter paint ppopulations.haps.radpainter
./finestructure -x 100000 -y 100000 -z 1000 populations.haps_chunks.out populations.haps_chunks.mcmc.xml
./finestructure -m T -x 10000 populations.haps_chunks.out populations.haps_chunks.mcmc.xml populations.haps_chunks.mcmcTree.xml

# for visualizing the results, "fineRADstructurePlot.R" was used; 
#fineRADstructurePlot.R is part of the fineRAD pipeline package we used without any modifications 

# References for fineRAD pipeline
# https://www.milan-malinsky.org/fineradstructure
# https://doi.org/10.1371/journal.pgen.1002453



