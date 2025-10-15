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
#cross_list=''
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
	#c) cross_list=${OPTARG};;
	#s) s_ref_path=${OPTARG};;
    #m) m_ref_path=${OPTARG};;
    #l) l_ref_path=${OPTARG};;
    esac
done


echo "get list of plate/well combinations:"
#only works becuase of current cross vs plate terminology
#the following line gets only plate well combos from experimental (coinfction) samples. the one below gets all
#grep "coinfection" ${sample_list} | awk -F"," '{print $1","$2}' > ${demuxed_path}/coinfection_plate_well.txt
grep "plate" ${sample_list} | awk -F"," '{print $1","$2}' > ${demuxed_path}/plate_well.txt
echo "plate well combos printed to plate_well.txt"
#grep "cross"  | awk -F"," '{print $1}' >> ${demuxed_path}/plaque_barcodes.txt



#prev scrip made a custom genotyping database for each cross, stored as "${demuxed_path}/cross/${cross}_database.fasta

#loop through each sample sequence data and u-search

for plaque in `cat ${demuxed_path}/plate_well.txt`;
do

	#get different variables from sample_list.csv
 	plate="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $1}')";
	well="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $2}')";
	coinfection="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $3}')";
	parent1="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $4}')";
    parent2="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $5}')";
	plaque_number="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $7}')";


echo
	echo "now on reverse primer ${plate}";
	echo "fwd primer ${well}";
	echo "which tags plaque ${plaque_number}";
	echo "from coinfection ${coinfection}";
	echo "with parents ${parent1} x ${parent2}";
	echo "and represents sample ${plaque_number}";
	
	#echo "reading file ${demuxed_path}/${plaque}/${plaque}.all.fastq ; generating file ${cross}/usearch/${cross}_${sample}_98_merged.b6"
	echo 

	#Lets organize the results by coinfection#/plate since one plate could contain plaques from multiple coinfections, or a coinfection could be split across plates
	echo "generating output folder ${coinfection}/${plate}"
	mkdir ${demuxed_path}/${coinfection}
	mkdir ${demuxed_path}/${coinfection}/${plate}

 	#Usearch files

	#instead of hard-coding seegment names, loop thru loci to usearch
	#can use this loop code to improve the rest of the script later :)
	for locus_fastq in ${demuxed_path}/${plate}/${well}/*.fastq;
	do
		locus_basepath="${locus_fastq##*/}";
		locus_slice1="${locus_basepath%%_cutadapt*}";
		locus="${locus_slice1##*_}";


		#USEARCH STRAIN ASSIGNMENT!!
		#Key change for this Oct 20 version is changing the database
		#usearch -usearch_global ${demuxed_path}/${plate}/${well}/${plate}_${well}_${locus}_cutadapt-lc_porechop.fastq -db ${demuxed_path}/${coinfection}/${coinfection}_parent_database_external.fasta -id 0.90 -blast6out ${demuxed_path}/${coinfection}/${plate}/${plate}_${well}_${locus}_90_merged.b6 -strand both

		#use this line instead to generate the results for comparing run 1 and 2 (uses sample name in output dfilename):
		usearch -usearch_global ${demuxed_path}/${plate}/${well}/${plate}_${well}_${locus}_cutadapt-lc_porechop.fastq -db ${demuxed_path}/${coinfection}/${coinfection}_parent_database_external.fasta -id 0.90 -blast6out ${demuxed_path}/${coinfection}/${plate}/${plate}_plaque${plaque_number}_${locus}_90_merged.b6 -strand both

	done

done


#output text editing and summary:
mkdir ${demuxed_path}/strain_assignment_output

#use first line for only coinfection samples, second for any sample type:
#grep "coinfection" ${sample_list} | awk -F"," '{print $3}' | sort | uniq > coinfection_list.txt
grep "plate" ${sample_list} | awk -F"," '{print $3}' | sort | uniq > experiment_list.txt


for experiment in `cat experiment_list.txt`;
do

#get list of plates, important so that the plate variable updates at the end when it spits out the data
#first gets only from data files. second jsut gets a list of all the plates from the sample plate file: 
ls -d ${demuxed_path}/${experiment}/plate* | sort | uniq | rev | cut -d '/' -f 1 | rev > plate_results.txt
#grep "plate" ${sample_list} | awk -F"," '{print $1"}' > plate_results.txt

	for plate in `cat plate_results.txt`;
	do
	#summarize the b6 results

        for b6_file in ${demuxed_path}/${experiment}/${plate}/plate*_*_*_90_merged.b6;
        do
		echo -n -e "${b6_file##*/}"\\t"";
        wc -l ${b6_file} | awk 'BEGIN { OFS = "\t"; ORS = "\t" } {if($1!="0") print $1}'
        rev ${b6_file} | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $11}' | rev | cut -d , -f 1 | sort | uniq -c | sort -nr | head -n 1 | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $1, $2} END {if (NR == 0) print ""}';
        done > ${demuxed_path}/strain_assignment_output/${experiment}_${plate}_strain_assignment_output.txt

	#end the plate loop
	done
	
#end the cross loop
done

echo "Amplicon processing and strain assignment done!"
