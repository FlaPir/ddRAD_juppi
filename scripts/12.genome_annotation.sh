############################
### redundancy reduction ### 
############################

#!/bin/bash

module load cd-hit/4.8.1

cd-hit -i proteins_for_braker.fasta \ #the input is the proteomes used in this study concatenated
       -o proteins_collapsed.fasta \ 
       -c 0.98 \ 
       -T 40 \ 
       -M 0 #the input is the proteomes used in this study concatenated

#########################
### genome annotation ### 
#########################

#!/bin/bash

source /etc/profile
source /etc/profile.d/modules.sh 2>/dev/null
module load singularity 2>/dev/null

THREADS=40
GENOME="./genome.fna"
ORIGINAL_PROT="./proteins_collapsed.fasta"
CLEAN_PROT="./proteins_cleaned.fasta"
SPECIES="Cx_pip"
WORKING_DIR="./braker_output"
SIF_IMAGE="/mnt/big/containers/singularity/braker3.sif"
MY_CONFIG="$HOME/augustus_config_custom"

mkdir -p $WORKING_DIR
rm -rf $MY_CONFIG/species/$SPECIES
sed 's/:/_/g' $ORIGINAL_PROT > $CLEAN_PROT #to replace ":" with "_" in the proteomes which can bother Augustus

export SINGULARITYENV_AUGUSTUS_CONFIG_PATH="$MY_CONFIG"

singularity exec -e \
    -B $PWD:$PWD,$MY_CONFIG:$MY_CONFIG \
    $SIF_IMAGE \
    braker.pl --genome=$GENOME \
              --prot_seq=$CLEAN_PROT \
              --species=$SPECIES \
              --threads=$THREADS \
              --workingdir=$WORKING_DIR \
              --softmasking \
              --gff3 \
              --overwrite \
              --AUGUSTUS_CONFIG_PATH=$MY_CONFIG

#BUSCO quality check (replace -i value with genome.fna for doing the genome quality check)
module load BUSCO/5.3.1 2>/dev/null

if [ -f "$WORKING_DIR/braker.aa" ]; then
    busco -i $WORKING_DIR/braker.aa -l diptera_odb10 -o busco_results -m prot -c $THREADS --force
fi

#############################
### functional annotation ### 
#############################

#!/bin/bash

#AGAT software to keep only the longest isoform (transcript) per gene 
agat_sp_keep_longest_isoform.pl -gff braker.gff3 \
                                -o genome_longest_isoform

#to extract the coding regions (CDS)
agat_sp_extract_sequences.pl -g genome_longest_isoform \
                             -f ../genome.fna \
                             -t cds \
                             -p \
                             -o CDS.faa

#DIAMOND for the functional annotation 
diamond makedb --in uniprotkb_proteome_Culex_pipiens_pipiens.fasta \
                -d db_culex_uniprot #uniprotkb_proteome_Culex_pipiens_pipiens.fasta -> UniProtKB: UP001562425

diamond blastp -q CDS.faa \ #coding regions only 
               -d db_culex_uniprot \
               -o Culex_pipiens_functional_annotation.tsv \ 
               --outfmt 6

#Selection of the best hits
sort -k1,1 -k11,11g Culex_pipiens_functional_annotation.tsv | sort -u -k1,1 > best_hits.tsv 

awk '$5 < 1e-10 && $3 > 30' best_hits.tsv > reliable_annotations.tsv #final file with only hits characterized by e-value < 1e-10 and %identity > 30 

#Mapping of the high-confidence functional features on the genome
agat_sp_manage_functional_annotation.pl -f genome.gff \
                                        -i reliable_annotations.tsv \
                                        --output Culex_pipiens_annotated_genome.gff3

