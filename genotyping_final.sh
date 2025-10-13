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
#s_ref_path=''
#m_ref_path=''
#l_ref_path=''


print_usage() {
  printf "Usage: ..."
}

while getopts d:e:c:s:m:l: flag
do
    case "${flag}" in
	d) demuxed_path=${OPTARG};;
	e) sample_list=${OPTARG};;
	c) cross_list=${OPTARG};;
	#s) s_ref_path=${OPTARG};;
    #m) m_ref_path=${OPTARG};;
    #l) l_ref_path=${OPTARG};;
    esac
done


#get list of barcodes for plaques only:
#only works becuase of current cross vs plate terminology
grep "cross" ${sample_list} | awk -F"," '{print $1}' >> ${demuxed_path}/plaque_barcodes.txt



#prev scrip made a custom genotyping database for each cross, stored as "${demuxed_path}/cross/${cross}_database.fasta

#loop through each sample sequence data and u-search


for plaque_barcode in `cat ${demuxed_path}/plaque_barcodes.txt`;
do

	#get different variables from sample_list.csv
 	plate="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $2}')"
	well="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $3}')"
	parent1="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $5}')";
        parent2="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $6}')"
	plaque_number="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $7}')"


echo
	echo "now on plaque ${plaque_barcode}, aka plate ${plate}, cross ${cross}, parents ${parent1} x ${parent2}; and plaque # ${plaque_number}"
	#echo "reading file ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq ; generating file ${cross}/usearch/${cross}_${sample}_98_merged.b6"
	echo 

	mkdir ${demuxed_path}/${cross}/${plate}

 	#Usearch files

	#instead of hard-coding seegment names, loop thru loci to usearch
	#can use this loop code to improve the rest of the script later :)
	for locus_fastq in ${demuxed_path}/${plaque_barcode}/*_reads.fastq;
	do
		locus_basepath="${locus_fastq##*/}";
		locus_slice1="${locus_basepath%%_reads*}";
		locus="${locus_slice1##*_}";


		#USEARCH STRAIN ASSIGNMENT!!
		#Key change for this Oct 20 version is changing the database
		usearch -usearch_global ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_${locus}_reads.fastq -db ${demuxed_path}/${cross}/${cross}_parent_database.fasta -id 0.90 -blast6out ${demuxed_path}/${cross}/${plate}/${cross}_${wellZ}_${locus}_90_merged.b6 -strand both

	done

done


#output text editing and summary:
mkdir ${demuxed_path}/strain_assignment_output

for cross in `cat ${cross_list}`;
do
        #summarize the b6 results

        for b6_file in ${demuxed_path}/${cross}/${plate}/cross*_*_*_90_merged.b6;
        do
		plate="${b6_file%%/cr*}"; 
  		plate="${plate##*/}"
        #echo "generating strain assignment output data in output folder for plate ${plate}, cross ${cross}!"
        echo -n -e "${b6_file##*/}"\\t"";
        wc -l ${b6_file} | awk 'BEGIN { OFS = "\t"; ORS = "\t" } {if($1!="0") print $1}'
        rev ${b6_file} | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $11}' | rev | cut -d , -f 1 | sort | uniq -c | sort -nr | head -n 1 | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $1, $2}'
        done > ${demuxed_path}/strain_assignment_output/${cross}_${plate}_strain_assignment_output.txt

#end the cross loop
done

echo "Amplicon processing and strain assignment done!"
