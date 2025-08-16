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
grep "parent" ${sample_list} | awk -F"," '{print $1}' > ${demuxed_path}/grid_barcodes.txt



#usearch in database with all potential single strains

#loop through each sample sequence data and u-search


for plaque_barcode in `cat ${demuxed_path}/grid_barcodes.txt`;
do

        #get different variables from sample_list.csv
        plate="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $2}')"
        #cross="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $3}')"
        parent1="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $4}')";
        #parent2="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $5}')"
        plaque_number="$(grep -m 1 ${plaque_barcode} ${sample_list} | awk -F"," '{print $6}')"


echo
        echo "now on plaque ${plaque_barcode}, aka plate ${plate}, strain ${parent1} and sample ${plaque_number} "
        #echo "reading file ${demuxed_path}/${plaque_barcode}/${plaque_barcode}.all.fastq ; generating file ${cross}/usearch/${cross}_${sample}_98_merged.b6"
        echo

        mkdir ${demuxed_path}/parent/${plate}

        #Usearch files

        #instead of hard-coding seegment names, loop thru loci to usearch
        #can use this loop code to improve the rest of the script later :)
        for locus_fastq in ${demuxed_path}/${plaque_barcode}/*_reads_filt.fastq;
        do
                locus_basepath="${locus_fastq##*/}";
                locus_slice1="${locus_basepath%%_reads*}";
                locus="${locus_slice1##*_}";

                echo "PLATE = ${plate}, barcode = ${plaque_barcode}"
                echo
                #USEARCH STRAIN ASSIGNMENT!!
                #Key change for this Oct 20 version is changing the database
                usearch -usearch_global ${demuxed_path}/${plaque_barcode}/${plaque_barcode}_${locus}_reads_filt.fastq -db ${demuxed_path}/parent/grid_database.fasta -id 0.99 -blast6out ${demuxed_path}/parent/${plate}/parent_${plate}_${plaque_barcode}_${plaque_number}_${locus}_90_merged.b6 -strand both
                echo
        done

done


#output text editing and summary:
mkdir ${demuxed_path}/parent_strain_assignment_output


for plate_name in `cat plate_list_grid.txt`; do

#       echo "entered final loop!";
#       #echo "CROSS = ${cross}";

        # echo "now entering b6_file for loop";
        #summarize the b6 results

        for b_file in ${demuxed_path}/parent/${plate_name}/parent_plate*_barcode*_*_90_merged.b6; do
#       echo "${b_file}";
        plate="${b_file%/parent*}";
        plate="${plate##*/}";
        #echo "${plate}";
       #echo "generating strain assignment output data in output folder for plate ${plate}, cross ${cross}!"
       #echo "the next code will output to the final results file, ${demuxed_path}/strain_assignment_output/${cross}_${plate}_strain_assignment_output.txt";
       #echo "b_file = ${b_file}";
        echo -n -e "${b_file##*/}"\\t"";
        wc -l ${b_file} | awk 'BEGIN { OFS = "\t"; ORS = "\t" } {if($1!="0") print $1}'
        rev ${b_file} | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $11}' | rev | cut -d , -f 1 | sort | uniq -c | sort -nr | head -n 1 | awk 'BEGIN { OFS = "\t"; ORS = "\n"} {print $1, $2}'
        done > ${demuxed_path}/parent_strain_assignment_output/parent_${plate_name}_strain_assignment_output.txt


        done

#end the cross loop
done

echo "Amplicon processing and strain assignment done!"
