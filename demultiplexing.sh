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
echo "Information on porechop trimming and read-splitting:"
grep "adapters" ${working_dir}/${reads_name}_porechop.log



#Ok lets make directories for the demuliplexing output:
#this makes a folder for every plate barcode and a subfolder for every well barcode. could later change to do based on the files output by cutadapt instead
for plate_dir in $(grep '^>' ${plate_barcodes} | sed 's/^>//'); 
do mkdir "${plate_dir}"; 
	for well_dir in $(grep '^>' ${well_barcodes} | sed 's/^>//'); 
	do mkdir "${plate_dir}/${well_dir}";


	

	
	done
done

#make directory for cutadapt outputs:
mkdir -p ${working_dir}/cutadapt_outputs
outputs="${working_dir}/cutadapt_outputs"

#DEMULTIPLEXING - STEP 1 - PLATE
#use plate_barcodes.fasta file to search and demultiplex PLATE barcodes with cutadapt
cutadapt -a file:${plate_barcodes} -O 8 --action=lowercase --revcomp -e 0.15 -o ${outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${working_dir}/${reads_name}_porechop.fastq > ${outputs}/plate_${reads_name}_cutadapt_porechop.log

#Use find 
find "${outputs}" -type f -name 'plate?_cutadapt.fastq' | while read -r plate_file_path; do
  
	#Move the plate_ demultiplexed files we just made into directories based off the file names:
 	plate_file_name="$(basename "$plate_file_path")"
  	plate="${plate_file_name%%_*}"  
  	plate_dir="${outputs}/${plate}";   
 	mkdir -p "${plate_dir}"
  	#mkdir -p "${outputs}/${plate}";
	mv "${plate_file_path}" "${plate_dir}/";
	#mv "${plate_file_path}" "${outputs}/${plate}/";
	echo "File sorted into directory by plate ID."


  	#Ok we want to execute the well-demultiplexing step once for each plate file. so include it in this loop:

	#DEMULTIPLEXING - STEP 2 - WELL
	cutadapt -g file:${well_barcodes} -O 8 --action=lowercase --revcomp -e 0.15 -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_file_path} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log;

	#move plate??_well??_ demultiplexed files to well folders:
	#plate_dir is updated for each value of plate
	# Make sure it's a directory
  	[ -d "$plate_dir" ] || continue

	#Loop through matching files inside the plate directory
	for plate_well_file in "$plate_dir"/plate??_well??_*.fastq; do
    # Extract the filename
    plate_well_file_name="$(basename "$plate_well_file")"

    # Extract the well ID (e.g., well01)
    well="${plate_well_file_name#*_}"
    well="${well_id%%_*}"

    # Create the well subdirectory
    plate_well_dir="$plate_dir/$well"
    mkdir -p "$plate_well_dir"

    # Move the file into the well subdirectory
    mv "$plate_well_file" "$plate_well_dir/"
  done


	


  
  done




#DEMULTIPLEXING - STEP 2 - WELL

#okay demultiplex by plaque, input = plate-demuxed files

for plate_dir in $(grep '^>' ${plate_barcodes} | sed 's/^>//'); 
	do for fastq in ${working_dir}/${plate_dir}/plate*.fastq; 
		do plate_temp="${fastq##*/}"; 
		plate="${plate_temp%%_*}";
		cutadapt -g file:${well_barcodes} -O 8 --action=lowercase --revcomp -e 0.15 -o ${working_dir}/${plate}/{name}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${fastq} > ${working_dir}/${plate}/${plate}_well_${reads_name}_cutadapt_porechop.log;
	done

#ended here 10/19

#get read counts for each file
	#for fastq in *plaque*-demux*; do cross_temp="${fastq##*cutadapt-}";cross="${cross_temp%%-demux*}";count=$( wc -l ${fastq} | awk '{print $1 / 4}'); echo "${cross},${count}" >> file_counts.csv; done

for well_dir in $(grep '^>' ${well_barcodes} | sed 's/^>//'); 
#demux on primer, input = cross-plaque files
	for fastq in ${working_dir}/${plate_dir}/plate*well*.fastq; 
		do plate_temp="${fastq##*/}"; 
		plate="${plate_temp%%_*}";
		echo "now demuxing ${fastq}";
		cutadapt -a small=CTTTCGTACAACCGAGTAGG...CTCCTGAAGTATCTCACGCC -a medium=CGCTACGGCGGTATTGTC...GCTCACCAAGTAAGGTGTAGTAT -a large=TCGATGTTCAACTACTACGC...GCGAGACTCGCTTTGC -O 10 --action=lowercase --revcomp -e 0.15 -o cutadapt-${crplq}-{name}-demux_porechop_land-pads_with_revcomp.fastq ${fastq} > current_cutadapt-${crplq}-segment-demux_landpads_RC_log.txt;
	done





done













