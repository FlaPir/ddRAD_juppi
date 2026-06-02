#!/bin/bash

bedtools intersect -a TOP_pca_outliers.bed -b st_top_outliers.bed -wa -wb > step1_top_intersec_outliers.bed

bedtools intersect -a step1_top_intersec_outliers.bed -b Culex_pipiens_annotated_genome.gff3 -wa -wb > top_intersec_outliers_annotated.txt

awk -v OFS='\t' '{print $1, $3, $4, $5, $10, $13, $19}' top_intersec_outliers_annotated.txt > top_intersec_outliers_annotated_clean.tsv

#The final table was then supplemented with the results obtained from eggNOG-mapper (https://galaxy-main.usegalaxy.org/) to add GO terms, KEGG pathways, etc.



