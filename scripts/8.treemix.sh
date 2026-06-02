#!/bin/bash

  plink \
    --bfile datasetA/B \
    --allow-extra-chr
    --within treemix_cluster #text file LOCALITY /t SAMPLE /t LOCALITY
    --freq
    --nonfounders
    --out freq_treemix
    
#############################################################
## conversion to a Treemix compatible file of plink output ##
#############################################################

python3 /path/to/plink2treemix.py /path/to/freq_per_treemix.frq.strat | gzip > treemix.gz

# plink2treemix.py downloaded by https://github.com/thomnelson/tools/blob/master/plink2treemix.py and modified 

#!/usr/bin/env python3
import sys
import pandas as pd

if len(sys.argv) != 2:
    sys.stderr.write("\nUsage: plink2treemix.py plink.frq.strat\n")
    sys.stderr.write(" (script will write uncompressed output to stdout)\n\n")
    sys.exit()

plink_freq_file = sys.argv[1]

freq_df = pd.read_csv(plink_freq_file, sep=r'\s+')

snp_IDs = freq_df["SNP"].unique()
pop_IDs = freq_df["CLST"].unique()
nsnps = len(snp_IDs)
npops = len(pop_IDs)

sys.stderr.write(f"SNPs in file: {nsnps}\n")
sys.stderr.write(f"Populations in input (n = {npops}):\n")
for pop in pop_IDs:
    sys.stderr.write(f"  {pop}\n")
    
print(" ".join(pop_IDs))

nsnps_out = 0
for snp, snp_df in freq_df.groupby("SNP"):
    if len(snp_df) != npops:
        continue

    nsnps_out += 1
    if nsnps_out % 100 == 0:
        sys.stderr.write(f"SNPs written: {nsnps_out}\r")

    counts = []
    for pop in pop_IDs:
        row = snp_df[snp_df["CLST"] == pop].iloc[0]
        nalleles = int(row["NCHROBS"])
        minor_count = int(row["MAC"])
        major_count = nalleles - minor_count
        counts.append(f"{major_count},{minor_count}")

    print(" ".join(counts))

sys.stderr.write(f"\nFinished! Total SNPs written: {nsnps_out}\n")


#############
## Treemix ##
#############

#!/bin/bash

for root in LN_J LN BO MO RA; do
  for i in {0..5}; do
    echo "Running: TreeMix m=$i rooted on $root"
    treemix -i treemix.gz -m $i -o treemix_${root}_m$i -root $root -k 500 -bootstrap > log_${root}_m$i.txt
  done
done


#unrooted cycle
for i in {0..5}; do
  echo "Running: TreeMix m=$i (Unrooted)"
  treemix -i treemix.gz -m $i -o treemix_unrooted_m$i -k 500 -bootstrap > log_unrooted_m$i.txt
done



