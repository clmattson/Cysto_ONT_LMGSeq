#!/bin/bash


#input to gather: 

#d path to all the barcode folders ('demuxed')
#e sample list  - CSV(!!) file (wih path) with all samples: ie barcode, cross, parent 1, parent 2
#c cross_list with path
#s S genotyping locus reference path
#m M genotyping locus reference path
#l L genotyping locus reference path


#fast5_pass_path=''
demuxed_path=''
sample_list=''
cross_list=''
s_ref_path=''
m_ref_path=''
l_ref_path=''


print_usage() {
  printf "Usage: ..."
}

while getopts d:e:c:s:m:l: flag
do
    case "${flag}" in
	d) demuxed_path=${OPTARG};;
	e) sample_list=${OPTARG};;
	c) cross_list=${OPTARG};;
	s) s_ref_path=${OPTARG};;
        m) m_ref_path=${OPTARG};;
        l) l_ref_path=${OPTARG};;
    esac
done

