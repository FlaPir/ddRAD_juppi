#!/bin/bash

#before running ADMIXTURE it's important to change the .bim file this way: awk '{ $1 = NR; print $0 }' old.bim > new.bim
#it will change the first column from the complete names of the chromosomes to numbers (OZ004311.1 --> 1)

INPUT_FILE="/path/to/file.bed"
BASE_NAME=$(basename "$INPUT_FILE" .bed)
OUTPUT_DIR=$(dirname "$INPUT_FILE")

cd "$OUTPUT_DIR" || exit

REPLICATES=20 #number of replicates per single K value

for K in {1..6}; do
    for R in $(seq 1 $REPLICATES); do

        #-s to change seed for each run
        path/to/the/software/admixture32 "$INPUT_FILE" $K --cv -s time | tee "${OUTPUT_DIR}/log_K${K}_run${R}.txt"

        #changes the name of the files to avoid overlapping with the next run 
        mv "${BASE_NAME}.${K}.Q" "${BASE_NAME}_K${K}_run${R}.Q"
        mv "${BASE_NAME}.${K}.P" "${BASE_NAME}_K${K}_run${R}.P"

    done
done


##########
## PONG ##
##########

#cross-validation error 
grep "CV error" log_K*.txt | sed 's/log_K//; s/_run/ /; s/.txt:CV error (K=//; s/):/ /' | awk '{print $1, $2, $3, $4}' > cv_data.txt

# https://github.com/ramachandran-lab/pong/tree/master

pong -m filemap \ #3 columns files K*_run* \t * \t K*_run*.Q
     -i ind_to_pop.txt \ #1 column file listing the population names repeated for the number of samples they represent  
     -n pop_order.txt \ #1 column file with the order in which you want to display the population in the graphic 
     -p 4000 #to specify the port of the server 



