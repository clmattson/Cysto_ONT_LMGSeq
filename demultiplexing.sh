#!/bin/bash


#input to gather: 

#d path to the desired working dir
#r path to reads  
#c plate_barcode.fasta file with path
#c well_barcode.fasta file with path
#s S genotyping locus reference path
#m M genotyping locus reference path
#l L genotyping locus reference path


demuxed_path=''
reads=''
plate_barcodes=''
well_barcodes=''
#s_ref_path=''
#m_ref_path=''
#l_ref_path=''


print_usage() {
  printf "Usage: ..."
}

while getopts d:r:p:w: flag
do
    case "${flag}" in
	d) working_dir=${OPTARG};;
	r) reads_path=${OPTARG};;
	p) plate_barcodes=${OPTARG};;
	w) well_barcodes=${OPTARG};;
	#s) s_ref_path=${OPTARG};;
  #m) m_ref_path=${OPTARG};;
  #l) l_ref_path=${OPTARG};;
    esac
done


#MOVE EXISTING OUTPUTS and RE-SET

#gets random number to move any old outputs to new temp folder:
random_number=$RANDOM

#move old outputs to random new folder so that they wont be over-written and old files wont interfere with new outputs:
echo "moving strain_assignment_output to strain_assignment_output_${random_number} and all *.b6 files in any plate*/well* folder to strain_assignment_output_${random_number}/b6_files"
mv ${demuxed_path}/strain_assignment_output ${demuxed_path}/strain_assignment_output_${random_number}
mkdir ${demuxed_path}/strain_assignment_output_${random_number}/b6_files
mv */*/*.b6 ${demuxed_path}/strain_assignment_output_${random_number}/b6_files



#Start analysis witg basecalled data from Dorado - run dorado with trimming DISABLED so that the adapters can be detected
#Reads should be in one large Fastq


We will use a combination of porechop and cutadapt to demulitplex the data: 

reads_name="${reads_path%.*}"
reads_name="${reads_path##*/}"

#PORECHOP - split reads on middle landing pads - also trims them
#porechop also needs to be run in the cutadapt conda environment
porechop -i ${reads_path} --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${working_dir}/${reads_name}_porechop.fastq > ${working_dir}/${reads_name}_porechop.log

#get porechop trimming information:

grep "adapters"  ${working_dir}/${reads_name}_porechop.log



#use plate_barcodes.fasta file to search and demultiplex PLATE barcodes with cutadapt

cutadapt -a file:${plate_barcodes} -O 8 --action=lowercase --revcomp -e 0.15 -o ${working_dir}/{name}_${reads_name}_porechop.fastq ${working_dir}/${reads_name}_porechop.fastq > ${working_dir}/{name}_${reads_name}_porechop.log


#use well_barcodes.fasta file to search and demultiplex PLATE barcodes with cutadapt

#okay demultiplex by plaque, input = cross files
for fastq in plate*.fastq; 
do cross_temp="${fastq##*cross-}";
cross="${cross_temp%%-demux*}";
cutadapt -g file:plaque_barcodes.fasta -O 8 --action=lowercase --revcomp -e 0.15 -o cutadapt-${cross}-{name}-demux_porechop_land-pads_with_revcomp.fastq ${fastq} > current_cutadapt-${cross}-plaque-demux_landpads_RC_log.txt;
done

#get read counts for each file
for fastq in *plaque*-demux*; do cross_temp="${fastq##*cutadapt-}";cross="${cross_temp%%-demux*}";count=$( wc -l ${fastq} | awk '{print $1 / 4}'); echo "${cross},${count}" >> file_counts.csv; done


#demux on primer, input = cross-plaque files
for fastq in *plaque*-demux*.fastq; 
do crplq_temp="${fastq##*cutadapt-}";
crplq="${crplq_temp%%-demux*}";
echo "now demuxing ${fastq}";
cutadapt -a small=CTTTCGTACAACCGAGTAGG...CTCCTGAAGTATCTCACGCC -a medium=CGCTACGGCGGTATTGTC...GCTCACCAAGTAAGGTGTAGTAT -a large=TCGATGTTCAACTACTACGC...GCGAGACTCGCTTTGC -O 10 --action=lowercase --revcomp -e 0.15 -o cutadapt-${crplq}-{name}-demux_porechop_land-pads_with_revcomp.fastq ${fastq} > current_cutadapt-${crplq}-segment-demux_landpads_RC_log.txt;
done













