#!/bin/bash



get_dir_timestamp() {
    local dir="$1"

    # get modify and change times
    local mtime=$(stat --format='%y' "$dir")
    local ctime=$(stat --format='%z' "$dir")

    # pick the earliest timestamp
    local earliest=$(printf "%s\n%s\n" "$mtime" "$ctime" | sort | head -n 1)

    # convert to safe filename format
    date -d "$earliest" +"%m-%d-%Y_%H.%M"
}


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
        l) lib_prep_database=${OPTARG};;
        #c) cross_list=${OPTARG};;
        #s) s_ref_path=${OPTARG};;
    #m) m_ref_path=${OPTARG};;
    #l) l_ref_path=${OPTARG};;
    esac
done

demuxed_path="${demuxed_path}/cutadapt_outputs"
samp_list_filename=$(basename "${sample_list}")
cp "${sample_list}" "${demuxed_path}/${samp_list_filename}"
sample_list="${demuxed_path}/${samp_list_filename}"

lib_prep_database_filename=$(basename "${lib_prep_database}")
cp "${lib_prep_database}" "${demuxed_path}/${lib_prep_database_filename}"
lib_prep_database="${demuxed_path}/${lib_prep_database_filename}"



# MOVE EXISTING OUTPUTS and RE-SET


###############################################
# RENAME strain_assignment_output
###############################################

shared_ts=""

# 1) If strain_assignment_output exists, use its timestamp
if [ -d "${demuxed_path}/strain_assignment_output" ]; then
    sa_dir="${demuxed_path}/strain_assignment_output"
    shared_ts=$(get_dir_timestamp "$sa_dir")

    echo "Renaming strain_assignment_output strain_assignment_output_${shared_ts}"
    mv "$sa_dir" "${demuxed_path}/strain_assignment_output_${shared_ts}"
    mkdir -p "${demuxed_path}/strain_assignment_output_${shared_ts}/b6_files"
fi


###############################################
# DETERMINE SHARED TIMESTAMP FOR EXPERIMENT DIRS
###############################################

# 2) If no timestamp yet, find earliest among experiment dirs
if [ -z "$shared_ts" ]; then
    earliest=""
    for d in coinfection positive misassigned negative; do
        if [ -d "${demuxed_path}/${d}" ]; then
            ts=$(get_dir_timestamp "${demuxed_path}/${d}")
            if [ -z "$earliest" ] || [[ "$ts" < "$earliest" ]]; then
                earliest="$ts"
            fi
        fi
    done
    shared_ts="$earliest"
fi

# If still empty, nothing exists do nothing
[ -z "$shared_ts" ] && return 0


###############################################
# RENAME EXPERIMENT DIRS USING shared_ts
###############################################

for d in coinfection positive misassigned negative; do
    if [ -d "${demuxed_path}/${d}" ]; then
        echo "Renaming ${d} to ${d}_${shared_ts}"
        mv "${demuxed_path}/${d}" "${demuxed_path}/${d}_${shared_ts}"
    fi
done


#output text editing and summary:
mkdir ${demuxed_path}/strain_assignment_output


reads_name=(${demuxed_path}/plate*.log)
reads_name="${reads_name%.*}";
reads_name="${reads_name##*/plate_}";
echo "reads_name variable = ${reads_name}";


echo "get list of plate/well combinations:"
#only works becuase of current cross vs plate terminology
#the following line gets only plate well combos from experimental (coinfction) samples.
#grep "coinfection" ${sample_list} | awk -F"," '{print $1","$2}' > ${demuxed_path}/coinfection_plate_well.txt
#this version gets plate,well combiations for ALL REAL samples:
#grep "plate" ${sample_list} | awk -F"," '{print $1","$2}' > ${demuxed_path}/plate_well.txt

#the following line gets ALL combinations of plate and well possible, good for checking assignemnt during demultiplexing:
#ls -d "${demuxed_path}/plate*/well*" | sed 's/\//,/g' | rev | awk -F',' '{print $1"," $2}' | rev > ${demuxed_path}/plate_well.txt

ls -d "$demuxed_path"/plate*/well* \
  | awk -F'/' '{print $(NF-1) "," $(NF)}' \
  > "$demuxed_path/plate_well.txt"

echo "we made it past plate_well.txt"



#prev scrip made a custom genotyping database for each cross, stored as "${demuxed_path}/cross/${cross}_database.fasta

#loop through each sample sequence data and u-search

for plaque in `cat ${demuxed_path}/plate_well.txt`;
do

        echo "we made it into the main for loop, currently on 'plaque' ${plaque}"
        #get different variables from sample_list.csv
        #plate="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $1}')";
        #well="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $2}')";

        #get well from the list loop instead of the file, so that samples with reads assigned to them will
        #always be analysed, even if they are not "real" plate/well combos included on the sample sheet
        plate="${plaque%%,*}";
        well="${plaque#*,}";

        echo "successfully assigned plate and well: ${plate} ; ${well}"
        #get other variables from sample_list.csv
        coinfection="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $3}')";
        parent1="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $4}')";
        parent2="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $5}')";
        plaque_number="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $7}')";

        echo "This sample was assigned type of: ${coinfection}; as per the value in the 3rd column of the sample sheet"


        #if coinfection is empty that means that the plate/well combination being analysed wasnt included in the sample sheet
        #the following line sets coinfection = "misassigned" if coinfection is empty.
        #if coinfection is not empty, it retains the existing value (set a few lines above)
        [ -z "$coinfection" ] && coinfection="misassigned"

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

        case "$coinfection" in
                positive)
                        #generate databse for positive controls/parents - will only contain the actual/correct parent fastas
                        #the meaningful datapoint for a usearch against this database will be the percent id, NOT the assigned identity or genotype
                        #database_file="${parent1}.fasta"
                        if [ -e "${demuxed_path}/${parent1}.fasta" ]; then
                                echo "Positive control reference for ${parent1} exists."
                        else
                                echo "making pos ctro ${parent1} database."
                                grep -A1 "${parent1}" ${lib_prep_database} | grep -v "^--$" > ${demuxed_path}/${parent1}.fasta

                        fi
                        database_file="${parent1}.fasta"


                        ;;
                negative)
                         #generate databse for positive controls/parents - will only contain the actual/correct parent fastas
                        #the meaningful datapoint for a usearch against this database will be the percent id, NOT the assigned identity or genotype
                        #database_file="${parent1}.fasta"

                        database_file="${lib_prep_database_filename}"

                        ;;

                misassigned)
                         #generate databse for positive controls/parents - will only contain the actual/correct parent fastas
                        #the meaningful datapoint for a usearch against this database will be the percent id, NOT the assigned identity or genotype
                        #database_file="${parent1}.fasta"

                        database_file="${lib_prep_database_filename}"


                        ;;
                coinfection)
                        if [ -e "${demuxed_path}/${parent1}_${parent2}.fasta" ]; then
                                echo "Coinfection reference database for ${parent1} times ${parent2} exists."
                        else
                                echo "making coinf ${parent1} times ${parent2}  database."
                                grep -A1 "${parent1}" ${lib_prep_database} | grep -v "^--$" > ${demuxed_path}/${parent1}_${parent2}.fasta
                                grep -A1 "${parent2}" ${lib_prep_database} | grep -v "^--$" >> ${demuxed_path}/${parent1}_${parent2}.fasta

                        fi
                        database_file="${parent1}_${parent2}.fasta"

                        ;;
                *)
                        echo "Unknown type: ${coinfection}"
                        ;;
        esac



        #Usearch files

        #instead of hard-coding seegment names, loop thru loci to usearch
        #can use this loop code to improve the rest of the script later :)
        #for locus_fastq in ${demuxed_path}/${plate}/${well}/*.fastq;
        for locus_fastq in ${demuxed_path}/${plate}/${well}/plate??_well??_${reads_name}.fastq;
        do
                locus_basepath="${locus_fastq##*/}";
                #toggle below line to search plate/well demuxed reads vs cysto segment demuxed reads:
                #locus=$(echo "$locus_basepath" | sed -E "s/^plate[0-9]{2}_well[0-9]{2}_([a-z]+)_${reads_name}\.fastq$/\1/");
                locus=$(echo "$locus_basepath" | sed -E "s/^plate[0-9]{2}_well[0-9]{2}_${reads_name}\.fastq$/\1/");
                #locus_slice1="${locus_basepath%%_cutadapt*}";
                #locus="${locus_slice1##*_}";

                echo "reads_name variable = ${reads_name}";

                #USEARCH STRAIN ASSIGNMENT!!
                #Key change for this Oct 20 version is changing the database
                #usearch -usearch_global ${demuxed_path}/${plate}/${well}/${plate}_${well}_${locus}_${reads_name}.fastq -db ${demuxed_path}/${coinfection}/${coinfection}_parent_database_external.fasta -id 0.90 -blast6out ${demuxed_path}/${coinfection}/${plate}/${plate}_${well}_${locus}_90_merged.b6 -strand both

                #THIS LINE CANNOT GENERATE DATA FOR NEGATIVE CONTROLS OR MISASSIGNED SAMPLES
                #use this line instead to generate the results for comparing run 1 and 2 (uses sample name in output dfilename):

                #OLD USEARCH LINE:
                #usearch -usearch_global ${demuxed_path}/${plate}/${well}/${plate}_${well}_${locus}_${reads_name}.fastq -db ${demuxed_path}/${coinfection}/${coinfection}_parent_database_external.fasta -id 0.90 -blast6out ${demuxed_path}/${coinfection}/${plate}/${plate}_plaque${plaque_number}_${locus}_90_merged.b6 -strand both

                echo
                echo "Running USEARCH on ${plate}_${well}_${reads_name}.fastq with DB: ${demuxed_path}/${database_file} and outputting to: ${demuxed_path}/${coinfection}/${plate}/${plate}_${well}_90_merged.b6 "
                echo
                usearch -usearch_global ${demuxed_path}/${plate}/${well}/${plate}_${well}_${reads_name}.fastq -db ${demuxed_path}/${database_file} -id 0.90 -blast6out ${demuxed_path}/${coinfection}/${plate}/${plate}_${well}_90_merged.b6 -strand both -top_hit_only
        done

        pars="${database_file%.fasta}"

        fastq="${demuxed_path}/${plate}/${well}/${plate}_${well}_${reads_name}.fastq"
        total_reads=$(($(wc -l < "$fastq") / 4))


        out="${demuxed_path}/strain_assignment_output/${experiment}_${plate}_strain_assignment_output.txt"

        base="${coinfection}_${pars}_${plate}"

        for b6_file in "${demuxed_path}/${coinfection}/${plate}/${plate}_${well}"_90_merged.b6;
        do
                awk -v file="${b6_file##*/}" \
                -v base="${base}" \
                -v outdir="${demuxed_path}/strain_assignment_output" \
                -v total_reads="${total_reads}" '
    {
        # Find the field that looks like CA68_M, PHI6_L, etc.
        seg_full = ""
        for (i = 1; i <= NF; i++) {
            if ($i ~ /^[A-Za-z0-9]+_[A-Za-z0-9]+$/) {
                seg_full = $i
                break
            }
        }

        # If we didn t find a segment, hard fail
        if (seg_full == "") {
            printf "FATAL: no segment-like field found in line:\n%s\n", $0 > "/dev/stderr"
            exit 1
        }

        split(seg_full, a, "_")
        segment   = a[2]
        match_val = seg_full

        total[segment]++
        count[segment, match_val]++

        # track how many hits each read has (query ID is $1)
        read_hits[$1]++
    }

    END {
        sum_max = 0

        # First pass: find max per segment and accumulate sum_max
        for (seg in total) {
            max_match = ""
            max_count = 0

            for (key in count) {
                split(key, p, SUBSEP)
                if (p[1] == seg && count[key] > max_count) {
                    max_count = count[key]
                    max_match = p[2]
                }
            }

            max_for_seg[seg]   = max_count
            match_for_seg[seg] = max_match
            sum_max           += max_count
        }

        # compute unsorted reads
        unsorted = total_reads - sum_max

        # compute multi-hit reads
        multi_hit = 0
        for (r in read_hits)
            if (read_hits[r] > 1)
                multi_hit++

        # second pass: print output for each segment
        outfile = outdir "/" base "_strain_assignment_output.txt"

        for (seg in total) {
            printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\n",
                   file, seg, match_for_seg[seg],
                   max_for_seg[seg], total[seg],
                   unsorted, multi_hit >> outfile
        }
    }
        ' "$b6_file"

    echo "DEBUG: base='$base' well='$well' b6_file='$b6_file'"
done

done
