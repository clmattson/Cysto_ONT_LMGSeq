#!/bin/bash

# Demultiplexing Script for Long Read LMGSeq
# Script: demultiplexing.sh
# This file is shell script that downloads and processes raw sequencing reads
# Requires three arguments that should be passed from the command line:
#    run information, cross information, and experiment folders

#This script is part of the following manuscript:
#
# Courtney Mattson | Hirva Shah | Ivy José | Samuel L. Díaz Muñoz

#Verify we have 3 arguments

if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments: run, plate list, and experiment"
    echo "run has format 'runA'"
    echo "File name for plate list. Cross list file is just a plain text file one line per plate: plate1 [break] plate9 etc. "
    echo "Experiment can be any folder name you desire, e.g. 'infection_conditions'"
    exit
fi

#Checking command line parameters:
echo "Checking command line arguments: run, plate list, and experiment"
echo $1
echo $2
echo $3

#Download FASTQ sequencing files from SRA via EBI
#wget ftp://ftp.sra.ebi.ac.uk/

#Unzip files
#gunzip data/runA_*.fastq.gz

#This script conducts demultiplexing for one run, which is indicated in the first command-line argument: runA,
# the second command line argument includes a plate list (in this case detailing the experimental coinfections or crosses),
# and the third command line argument groups results into a directory (in this case we call this script twice to separate control experiments and pairwise experimental coinfections)

#Change into working directory, which in this case is the base directory where the script resides
BASEDIR=$(dirname "$0")
cd $BASEDIR

#Make a directory to hold all the intermediate files etc.
INTERMED_DIR="$1""_""$3"
mkdir $INTERMED_DIR

#Let's copy all of the shared directory
cp -r $BASEDIR/shared $INTERMED_DIR/

## First Copy Barcodes only for relevant crosses
cat $2 | grep -f - $BASEDIR/shared/plate_barcodes.fasta -A1 -w --no-group-separator > $INTERMED_DIR/plate_barcodes_selected.fasta
#Repeat with reverse barcodes
cat $2 | grep -f - $BASEDIR/shared/plate_barcodes_rev.fasta -A1 -w --no-group-separator > $INTERMED_DIR/plate_barcodes_selected_rev.fasta

#Now, copy entire files for 3' barcodes because we want all samples 
cp $BASEDIR/shared/well_barcodes.fasta $INTERMED_DIR/well_barcodes.fasta
cp $BASEDIR/shared/well_barcodes_rev.fasta $INTERMED_DIR/well_barcodes_rev.fasta
#Might change script to allow subset of wells, but unsure how to do combinatorially yet
#cat $2 | grep -f - $BASEDIR/shared/well_barcodes.fasta -A1 -w --no-group-separator > $INTERMED_DIR/well_barcodes_selected.fasta
#Repeat with reverse barcodes
#cat $2 | grep -f - $BASEDIR/shared/well_barcodes_rev.fasta -A1 -w --no-group-separator > $INTERMED_DIR/shared/well_barcodes_rev.fasta



#Change into Run directory
cd $INTERMED_DIR

#Load libraries and modules
module load conda/base
conda activate lmgseq_longread
module load usearch/11.0.667


#Remove adapters from all reads
#cutadapt -g TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT -g CCTGTACTTCGTTCAGTTACGTATTGC -e 0.4 -o trim_adapter.fastq /home/sldmunoz/LMGSeq_longread/data/big_subset.fastq --cores=0 --untrimmed-output=untrimmed_adapter.fastq > trim_adapter_log.txt
cutadapt -g TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT -g CCTGTACTTCGTTCAGTTACGTATTGC -e 0.4 -o trim_adapter.fastq /group/sldmunozgrp/data/cysto_LMGSeq_test_082025/FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq --cores=0 --untrimmed-output=untrimmed_adapter.fastq > trim_adapter_log.txt
grep "Reads with adapters:" trim_adapter_log.txt

#Remove landing pads from all reads
cutadapt -a ^TTGGTGCTGATATTGC...GAAGATAGAGCGACAG$ -a ^CTGTCGCTCTATCTTC...GCAATATCAGCACCAA$ -e 0.6 -o trim_landing_pad.fastq trim_adapter.fastq --cores=0 --untrimmed-output=untrimmed_landing_pad.fastq > trim_landing_pad_log.txt
grep "Reads with adapters:" trim_landing_pad_log.txt

#Now need to demultiplex, wells first
cutadapt -a file$:well_barcodes_rev.fasta -g ^file:well_barcodes.fasta -e 0 -o {name}.fastq trim_landing_pad.fastq --cores=0 --untrimmed-output={name}_untrimmed.fastq > well_log.txt
grep "Reads with adapters:" well_log.txt

#Now separate the wells
for file in well*.fastq;
do
    well=${file%.fastq} #Get file prefix
    echo "$file"
    echo "$well"

#Now demultiplexing and trimming the samples
cutadapt -a file$:plate_barcodes_selected.fasta -g ^file:plate_barcodes_selected_rev.fasta -e 0 -o {name}_${well}.fastq ${well}.fastq --cores=0 --untrimmed-output=${well}_untrimmed.fastq > ${well}_log.txt
grep "Reads with adapters:" ${well}_log.txt
done

echo "Demultiplexing done!"
