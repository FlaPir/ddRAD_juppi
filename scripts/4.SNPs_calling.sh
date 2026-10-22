#!/bin/bash

#Stacks population filtering for Dataset A
 
  populations \ #populations is a module in Stacks package
    -P /path/to/gstacks_output \
    -M /path/to/pop_map.txt \ #file with sample_id \t population_locality format 
    -O /path/to/output_folder \
    -t 16 #chosen number of cores 
    --min-gt-depth 10 \
    -R 0.6 \
    --min-maf 0.05 \
    --max-obs-het 0.7 \
    --vcf \
    --radpainter \
    --treemix \
    --plink 

#Stacks population filtering for Dataset C

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
