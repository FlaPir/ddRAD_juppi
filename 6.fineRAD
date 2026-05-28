#!/bin/bash

RADpainter hapsFromVCF datasetC.vcf > per_radpainter.haps.radpainter
RADpainter paint per_radpainter.haps.radpainter
finestructure -x 100000 -y 100000 -z 1000 per_radpainter.haps_chunks.out per_rad_painter_chunks.mcmc.xml
finestructure -m T -x 10000 per_radpainter.haps_chunks.out per_rad_painter_chunks.mcmc.xml per_radpainter.haps_chunks.mcmcTree.xml

# References for fineRAD pipeline
# https://www.milan-malinsky.org/fineradstructure
# https://doi.org/10.1371/journal.pgen.1002453


# for visualizing the results, "fineRADstructurePlot.R" was used; fineRADstructurePlot.R is part of the fineRAD pipeline package 
