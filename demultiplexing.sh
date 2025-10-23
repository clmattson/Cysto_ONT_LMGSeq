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

module load conda
conda activate cutadapt


#MOVE EXISTING OUTPUTS and RE-SET

#gets random number to move any old outputs to new temp folder:
random_number=$RANDOM

#move old outputs to random new folder so that they wont be over-written and old files wont interfere with new outputs:
echo "moving cutadapt_outputs and porechop_outputs directory to cutadapt_outputs_${random_number} & porechop_outputs_${random_number} and and all files in any plate*/well* folder"
#wont do anything if these havent been made yet

[ -d "${working_dir}/cutadapt_outputs" ] && mv ${working_dir}/cutadapt_outputs ${working_dir}/cutadapt_outputs_${random_number}
[ -d "${working_dir}/porechop_outputs" ] && mv ${working_dir}/porechop_outputs ${working_dir}/porechop_outputs_${random_number}


#Start analysis witg basecalled data from Dorado - run dorado with trimming DISABLED so that the adapters can be detected
#Reads should be in one large Fastq


#We will use a combination of porechop and cutadapt to demulitplex the data: 

#make directory for cutadapt outputs:
mkdir -p ${working_dir}/porechop_outputs
porechop_outputs="${working_dir}/porechop_outputs"

echo "made directoy for porechop outputs"

reads_name="${reads_path%.*}";
reads_name="${reads_name##*/}";
echo "reads_name variable = ${reads_name}";

#PORECHOP - split reads on middle landing pads - also trims them
#porechop also needs to be run in the cutadapt conda environment
#it also seems like porechop seems to have memory issues when run on large files, lets add some parallel functionality:

#number of parts to split the fastq into
num_fastq_sections=5
#total number of reads in the FASTQ file
fastq_total_reads=$(wc -l ${reads_path} | awk '{print $1 / 4}')

#reads (sets of 4 lines) per part
lines_per_section=$(( (fastq_total_reads / num_fastq_sections) * 4 ))

#Use split to divide the file
split -l "${lines_per_section}" -d --additional-suffix=.fastq "${reads_path}" "${working_dir}/${reads_name}_chunk"


#porechop -i ${reads_path} --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${porechop_outputs}/${reads_name}_porechop.fastq > ${porechop_outputs}/${reads_name}_porechop.log


porechop -i ${working_dir}/${reads_name}_chunk00.fastq --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${porechop_outputs}/${reads_name}_porechop_chunk00.fastq > ${porechop_outputs}/${reads_name}_porechop_chunk00.log  
porechop -i ${working_dir}/${reads_name}_chunk01.fastq --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${porechop_outputs}/${reads_name}_porechop_chunk01.fastq > ${porechop_outputs}/${reads_name}_porechop_chunk01.log 
porechop -i ${working_dir}/${reads_name}_chunk02.fastq --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${porechop_outputs}/${reads_name}_porechop_chunk02.fastq > ${porechop_outputs}/${reads_name}_porechop_chunk02.log 
porechop -i ${working_dir}/${reads_name}_chunk03.fastq --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${porechop_outputs}/${reads_name}_porechop_chunk03.fastq > ${porechop_outputs}/${reads_name}_porechop_chunk03.log 
#porechop -i ${working_dir}/${reads_name}_chunk04.fastq --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${porechop_outputs}/${reads_name}_porechop_chunk04.fastq > ${porechop_outputs}/${reads_name}_porechop_chunk04.log
#porechop -i ${working_dir}/${reads_name}_chunk05.fastq --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${porechop_outputs}/${reads_name}_porechop_chunk05.fastq > ${porechop_outputs}/${reads_name}_porechop_chunk05.log

wait
echo "All porechop chunks complete!" 
echo "Now pasting all porechopped chunks back into one file: ${porechop_outputs}/${reads_name}_porechop.fastq"
cat ${porechop_outputs}/${reads_name}_porechop*chunk??.fastq >> ${porechop_outputs}/${reads_name}_porechop.fastq

echo "executed porechop for splitting reads on landing pads"

#get porechop trimming information:
echo
echo "Information on porechop trimming and read-splitting" for each section: 
grep "adapters" ${porechop_outputs}/${reads_name}_porechop_chunk*.log
echo



#Ok lets make directories for the demuliplexing output:
#this makes a folder for every plate barcode and a subfolder for every well barcode. could later change to do based on the files output by cutadapt instead
#for plate_dir in $(grep '^>' ${plate_barcodes} | sed 's/^>//'); 
#do mkdir "${plate_dir}"; 
#	for well_dir in $(grep '^>' ${well_barcodes} | sed 's/^>//'); 
#	do mkdir "${plate_dir}/${well_dir}";	
	#done
#done

#make directory for cutadapt outputs:
mkdir -p ${working_dir}/cutadapt_outputs
cutadapt_outputs="${working_dir}/cutadapt_outputs"

#DEMULTIPLEXING - STEP 2 - PLATE
#use plate_barcodes.fasta file to search and demultiplex PLATE barcodes with cutadapt. Higher -O
cutadapt -a file:${plate_barcodes} -O 14 --revcomp -e 0.15 -o ${cutadapt_outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${porechop_outputs}/${reads_name}_porechop.fastq > ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop.log
#cutadapt -a file:${plate_barcodes} -O 14 --action=lowercase --revcomp -e 0.15 -o ${cutadapt_outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${working_dir}/${reads_name}_porechop.fastq > ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop.log

echo "completed demultiplexing step 2 - cutadapt plate identification!"

#Use find 
find "${cutadapt_outputs}" -type f -name 'plate??_*.fastq' | while read -r plate_file_path; do
	echo
	echo "entered cutadapt loop, plate_file_path = ${plate_file_path}"

	#Move the plate_ demultiplexed files we just made into directories based off the file names:
 	plate_file_name="$(basename "$plate_file_path")"
  	plate="${plate_file_name%%_*}"  
  	plate_dir="${cutadapt_outputs}/${plate}";   
 	mkdir -p "${plate_dir}"
  	#mkdir -p "${cutadapt_outputs}/${plate}";
	mv "${plate_file_path}" "${plate_dir}/";
	#mv "${plate_file_path}" "${cutadapt_outputs}/${plate}/";
	echo "File sorted into directory 'plate_dir' : ${plate_dir}   "
	echo

  	#Ok we want to execute the well-demultiplexing step once for each plate file. so include it in this loop:

	#DEMULTIPLEXING - STEP 3 - WELL
	echo
	echo "executing demultiplaexing step 3: cutadapt search for well barcodes! input file=${plate_dir}/${plate_file_name}";
	cutadapt -g file:${well_barcodes} -O 14 --revcomp -e 0.15 -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_dir}/${plate_file_name} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log;
	#cutadapt -g file:${well_barcodes} -O 14 --action=lowercase --revcomp -e 0.15 -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_file_path} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log;

	#move plate??_well??_ demultiplexed files to well folders:
	#plate_dir is updated for each value of plate
	# Make sure it's a directory
  	[ -d "$plate_dir" ] || continue

	#Loop through matching files inside the plate directory
	for plate_well_file_path in "${plate_dir}"/plate??_well??_*.fastq; do
    	# Extract the filename
    	plate_well_file_name="$(basename "$plate_well_file_path")";

    	# Extract the well ID (e.g., well01)
    	well="${plate_well_file_name#*_}";
    	well="${well%%_*}";
		echo "current well = ${well}";

    	# Create the well subdirectory
    	plate_well_dir="$plate_dir/$well";
    	mkdir -p "$plate_well_dir";

    	# Move the file into the well subdirectory
		echo "moving ${plate_well_file_path} to ${plate_well_dir}/"
    	mv "${plate_well_file_path}" "${plate_well_dir}/";

		#while still looping thru values of plate and well, do cutadapt search for primers
		
		#DEMULTIPLEXING - STEP 4 - SEGMENT
		#okay demultiplex by plaque, input = plate-demuxed files; -O is smaller bc the primers are shorter
		cutadapt -g file:${well_barcodes} -O 10 --revcomp -e 0.15 -o ${plate_well_dir}/${plate}_${well}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_well_dir}/${plate_well_file_name} > ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop.log;
		#cutadapt -g file:${well_barcodes} -O 10 --action=lowercase --revcomp -e 0.15 -o ${plate_well_dir}/${plate}_${well}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_well_file_path} > ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop.log;
	
		#get read count in each file:
		for fastq in ${plate_well_dir}/${plate}_${well}_*_${reads_name}_cutadapt_porechop.fastq;
			do count=$( wc -l ${fastq} | awk '{print $1 / 4}'); 
			echo "${cross},${count}" >> ${working_dir}/file_counts.csv;
		done	
 	done
done

echo "Dont forget!! I moved your previous cutadapt_outputs and porechop_outputs directories to cutadapt_outputs_${random_number} & porechop_outputs_${random_number} and and all files in any plate*/well* folder :)"



