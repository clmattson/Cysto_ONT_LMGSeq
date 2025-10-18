#!/bin/bash


#input to gather: 

#d path to the data and working dir
#e sample list  - CSV(!!) file (wih path) with all samples: ie barcode, cross, parent 1, parent 2
#c cross_list with path
#s S genotyping locus reference path
#m M genotyping locus reference path
#l L genotyping locus reference path


#fast5_pass_path=''
demuxed_path=''
sample_list=''
#cross_list=''
#s_ref_path=''
#m_ref_path=''
#l_ref_path=''


print_usage() {
  printf "Usage: ..."
}

while getopts d:e: flag
do
    case "${flag}" in
	d) data_path=${OPTARG};;
	e) sample_list=${OPTARG};;
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







