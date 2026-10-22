
########################################
# 1.process_radtags                    #
# Author: Flavio Piras                 #
# Reference: ....                      #
#                                      #
########################################


#!/bin/bash

INPUT_DIR="/path/to/raw_reads"
OUTPUT_DIR="/path/to/process_radtags_output"

mkdir -p "$OUTPUT_DIR"

export INPUT_DIR OUTPUT_DIR

run_stacks() {
    R1="$1"
    SAMPLE=$(basename "$R1" .1.fq)
    R2="$INPUT_DIR/${SAMPLE}.2.fq"

    if [[ -f "$R2" ]]; then
        
        mkdir -p "$OUTPUT_DIR/$SAMPLE"

        process_radtags -1 "$R1" -2 "$R2" -o "$OUTPUT_DIR/$SAMPLE" \ #process_radtags is a module in the Stacks pipeline 
        --renz_1 sphI --renz_2 dpnII -c -q -r --retain_header --paired

        mv "$OUTPUT_DIR/$SAMPLE"/process_radtags.log "$OUTPUT_DIR/${SAMPLE}_process_radtags.log"
        mv "$OUTPUT_DIR/$SAMPLE"/*.fq* "$OUTPUT_DIR/" 2>/dev/null
        
        # Rimuoviamo la cartellina temporanea vuota
        rmdir "$OUTPUT_DIR/$SAMPLE"

}

export -f run_stacks
# to process 8 samples per run
parallel -j 8 run_stacks ::: "$INPUT_DIR"/*.1.fq
