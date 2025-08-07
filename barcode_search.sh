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

total_reads=$( wc -l FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq | awk '{print $1/4}')

for plaque_barcode_seq in $(grep -v "^>" "$plaque_barcodes_fasta" ); do
  seq_line=$(sed -n "/${plaque_barcode_seq}/=" ${plaque_barcodes_fasta});   
  echo "seq_line = ${seq_line}";
  header_line=$((seq_line -1));   
  reads_with_barcode=$(grep -c "${plaque_barcode_seq}" ${reads_fastq});
  echo "reads_with_barcode = ${reads_with_barcode}";
  #total_reads=$( wc -l FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq | awk '{print $1/4}');
  echo "total_reads = ${total_reads}";
  percent_reads_with_barcode=$(echo "scale=5 ; $reads_with_barcode/$total_reads"| bc);
  header=$(sed -n "${header_line}p" ${plaque_barcodes_fasta}); 
  echo "${percent_reads_with_barcode}% of the reads contain an exact match to the barcode for ${header}";
  echo;
  echo;
done




  






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



