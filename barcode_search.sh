#!/bin/bash

#to do:


#get input:

#

#-p ${parent_dir} directory containing data and (#sup_basecall)
#-d sample_suffix - starts with an underscore, example: "_norm3_meta_spades"
#-s ${host_s} S segment reference file & path
#-m ${host_m} M segment reference file & path
#-l ${host_l} L segment reference file & path


#parent_path=''
sample_suffix=''
#host_s=''
#host_m=''
#host_l=''


print_usage() {
  printf "Usage: ..."
}

while getopts d: flag
do
    case "${flag}" in
d) sample_suffix=${OPTARG};;
#s) host_s=${OPTARG};;
#m) host_m=${OPTARG};;
#l) host_l=${OPTARG};;

    esac
done


#cd ${parent_path} || exit

echo "working shell directory:"
pwd
echo




#okay % classified looking sad for the plaque primers. lets get some code that searches for the barcodes and identifies the % of reads containing an exact match

file="dorado_barcode_arrs/plaque_barcodes.fasta"
while read line; do
  echo "${line}"
done < "${file}


for plaque_barcode_seq in $(grep -v "^>" "$plaque_barcodes_fasta" ); do
  grep -c "${plaque_barcode_seq}" ${reads_fastq}






file="dorado_barcode_arrs/plaque_barcodes.fasta"
grep "^>" "$FASTA_FILE" | while IFS= read -r header; do
    # Remove the '>' character from the header
    sequence_id="${header#>"}"

    # You can now use $sequence_id to process each sequence
    echo "Processing sequence: $sequence_id"

done < "${$FASTA_FILE}"

    # Example: Extract the sequence data for the current ID (more complex for multi-line sequences)
    # This example assumes single-line sequences or you need more advanced parsing for multi-line
    # For robust parsing of multi-line sequences, consider tools like 'seqtk subseq' or 'samtools faidx'
    # or using a scripting language like Python with Biopython.
    # grep -A 1 "$header" "$FASTA_FILE" | tail -n 1
done



