#!/bin/bash


#input to gather: 

#d path to the desired working dir
#r path to reads  
#c cross_list with path
#s S genotyping locus reference path
#m M genotyping locus reference path
#l L genotyping locus reference path


demuxed_path=''
reads=''
#cross_list=''
#s_ref_path=''
#m_ref_path=''
#l_ref_path=''


print_usage() {
  printf "Usage: ..."
}

while getopts d:r: flag
do
    case "${flag}" in
	d) working_dir=${OPTARG};;
	r) reads_path=${OPTARG};;
	#c) cross_list=${OPTARG};;
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

#PORECHOP - split reads on middle adapters
porechop -i ${reads_path} --verbosity 2 --end_threshold 70 --middle_threshold 80 --extra_end_trim 0 --end_size 150 --min_split_read_size 200 --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 --min_trim_size 8 -o ${working_dir}/${reads_name}_porechop.fastq > ${working_dir}/${reads_name}_porechop.log












