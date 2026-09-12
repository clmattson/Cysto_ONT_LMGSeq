#!/bin/bash

#example command
#bash genotyping.sh -d /group/sldmunozgrp/cysto_LMGSeq08-25/sam_basecalled_flu_2025/CA09xPAN99_102325 -s /group/sldmunozgrp/cysto_LMGSeq08-25/sam_basecalled_flu_2025/CA09xPAN99_102325/sample_list_CA09xPAN99.csv -l /group/sldmunozgrp/cysto_LMGSeq08-25/sam_basecalled_flu_2025/CA09xPAN99_102325/all_flu_refs.fasta -b 0 2>&1 | tee /group/sldmunozgrp/cysto_LMGSeq08-25/sam_basecalled_flu_2025/CA09xPAN99_102325/genotyping_log_CA09xPAN99.txt
#(note: this script also writes its own log inside genotyping_outputs/ now, so the "| tee" above is optional, not required)

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

# quick fail-fast check so a forgotten module load shows up immediately with a
# clear message, instead of a confusing usearch-not-found error several steps in
for tool in usearch; do
    command -v "$tool" >/dev/null || { echo "Oops! '$tool' not found on PATH - did you 'module load usearch' before running this script?"; exit 1; }
done


echo "flag info - input to gather from user"

echo " d - path to the working directory that CONTAINS cutadapt_outputs - this is the SAME working directory you passed to "
echo " the demultiplexing script, NOT cutadapt_outputs itself. All genotyping outputs (coinfection/positive/negative/misassigned "
echo " .b6 files, strain_assignment_output summaries, and the parent reference databases) get written to a genotyping_outputs"
echo " directory here, alongside cutadapt_outputs - nothing gets written into cutadapt_outputs itself."

echo " s - sample list  - CSV(!!) file (wih path) with all samples: ie barcode, cross, parent 1, parent 2"
echo " l - library prep database = a fasta file with ALL possible reference genomes input into the entire library prep."
echo " Indiv sample datbases get built from this. Include one reference per strain per genome segment and name them like this:"
echo " >Strain_Seg - ex: >PAN99_M; >CA09_NP etc"

echo "b - do_blast: 0 or 1 only - enter 1 if you Want to blast the reads that usearch couldnt match to any ref. warning - its a bit slow!"

#fast5_pass_path=''
working_dir=''
sample_list=''
lib_prep_database=''
do_blast=''


print_usage() {
  printf "Usage: ..."
}

while getopts d:s:l:b: flag
do
    case "${flag}" in
        d) working_dir=${OPTARG};;
        s) sample_list=${OPTARG};;
        l) lib_prep_database=${OPTARG};;
        b) do_blast=${OPTARG};;
    esac
done

cutadapt_outputs="${working_dir}/cutadapt_outputs"
genotyping_outputs="${working_dir}/genotyping_outputs"


###############################################
# MOVE EXISTING genotyping_outputs (if any) to a timestamped backup, same
# pattern the demultiplexing script uses for cutadapt_outputs/porechop_outputs -
# one directory, one move, instead of hunting per-subfolder timestamps.
###############################################

if [ -d "${genotyping_outputs}" ]; then
    ts_geno=$(get_dir_timestamp "${genotyping_outputs}")
    echo "a genotyping_outputs directory already exists - moving existing genotyping_outputs directory to genotyping_outputs_from_${ts_geno}"
    mv "${genotyping_outputs}" "${genotyping_outputs}_from_${ts_geno}"
fi

mkdir -p "${genotyping_outputs}"
mkdir -p "${genotyping_outputs}/strain_assignment_output"

# this script's own log, so we can automatically save the output
genotyping_log="${genotyping_outputs}/genotyping_screen.log"
exec > >(tee -a "$genotyping_log") 2>&1
echo "Full genotyping-stage script output is being saved to: ${genotyping_log}"
echo "(usearch itself also writes a .b6 file per well inside genotyping_outputs/<sample_type>/<plate>/ - this log is the script's own narration, not a replacement for those)"

samp_list_filename=$(basename "${sample_list}")
cp "${sample_list}" "${genotyping_outputs}/${samp_list_filename}"
sample_list="${genotyping_outputs}/${samp_list_filename}"

lib_prep_database_filename=$(basename "${lib_prep_database}")
cp "${lib_prep_database}" "${genotyping_outputs}/${lib_prep_database_filename}"
lib_prep_database="${genotyping_outputs}/${lib_prep_database_filename}"


#output text editing and summary:

reads_name=(${cutadapt_outputs}/plate*.log)
reads_name="${reads_name%.*}";
reads_name="${reads_name##*/plate_}";
echo "reads_name variable = ${reads_name}";


echo "get list of plate/well combinations:"
#only works becuase of current cross vs plate terminology
#the following line gets only plate well combos from experimental (coinfction) samples.
#grep "coinfection" ${sample_list} | awk -F"," '{print $1","$2}' > ${genotyping_outputs}/coinfection_plate_well.txt
#this version gets plate,well combiations for ALL REAL samples:
#grep "plate" ${sample_list} | awk -F"," '{print $1","$2}' > ${genotyping_outputs}/plate_well.txt

#the following line gets ALL combinations of plate and well possible, good for checking assignemnt during demultiplexing:
#ls -d "${cutadapt_outputs}/plate*/well*" | sed 's/\//,/g' | rev | awk -F',' '{print $1"," $2}' | rev > ${genotyping_outputs}/plate_well.txt

ls -d "$cutadapt_outputs"/plate*/well* \
  | awk -F'/' '{print $(NF-1) "," $(NF)}' \
  > "$genotyping_outputs/plate_well.txt"

echo "we made it past plate_well.txt"



#prev scrip made a custom genotyping database for each cross, stored as "${genotyping_outputs}/cross/${cross}_database.fasta

#loop through each sample sequence data and u-search

for plaque in `cat ${genotyping_outputs}/plate_well.txt`;
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
    #sample_type describes what kind of sample this plate/well is (positive control, negative
    #control, coinfection, or misassigned) - renamed from "coinfection" since that name was
    #confusingly also one of the possible values it can hold
    sample_type="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $3}')";
    parent1="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $4}')";
    parent2="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $5}')";
    #added parent3/parent4 for quad/triple sample sheet format if needed
    parent3="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $6}')";
    parent4="$(grep -m 1 ${plaque} ${sample_list} | awk -F"," '{print $7}')";

    echo "This sample was assigned type of: ${sample_type}; as per the value in the 3rd column of the sample sheet"


    #if sample_type is empty that means that the plate/well combination being analysed wasnt included in the sample sheet
    #the following line sets sample_type = "misassigned" if sample_type is empty.
    #if sample_type is not empty, it retains the existing value (set a few lines above)
    [ -z "$sample_type" ] && sample_type="misassigned"

    echo
    echo "now on reverse primer ${plate}";
    echo "fwd primer ${well}";
    echo "which tags plaque ${plaque_number}";
    echo "from sample_type ${sample_type}";
    echo "with parents ${parent1} x ${parent2}${parent3:+ x ${parent3}}${parent4:+ x ${parent4}}";
    echo "and represents sample ${plaque_number}";

    #echo "reading file ${cutadapt_outputs}/${plaque}/${plaque}.all.fastq ; generating file ${cross}/usearch/${cross}_${sample}_98_merged.b6"
    echo

    #Lets organize the results by sample_type/plate since one plate could contain plaques from multiple sample types, or a sample type could be split across plates
    echo "generating output folder ${sample_type}/${plate}"
    mkdir -p ${genotyping_outputs}/${sample_type}
    mkdir -p ${genotyping_outputs}/${sample_type}/${plate}

    #check the value of sample_type and generate appropriate db for usearch
    case "$sample_type" in
            positive)
                    #generate databse for positive controls/parents - will only contain the actual/correct parent fastas
                    #the meaningful datapoint for a usearch against this database will be the percent id, NOT the assigned identity or genotype
                    #database_file="${parent1}.fasta"
                    if [ -e "${genotyping_outputs}/${parent1}.fasta" ]; then
                            echo "Positive control reference for ${parent1} exists."
                    else
                            echo "making pos ctro ${parent1} database."
                            grep -A1 "${parent1}" ${lib_prep_database} | grep -v "^--$" > ${genotyping_outputs}/${parent1}.fasta

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
                    #build database depending on how many parents are present (2, 3, or 4)
                    coinf_filename="${parent1}_${parent2}"

                    # [-n var] checks if the variable has value of non zero length, and then && executed the command that follows if the check before was true
                    [ -n "$parent3" ] && coinf_filename="${coinf_filename}_${parent3}"
                    [ -n "$parent4" ] && coinf_filename="${coinf_filename}_${parent4}"

                    if [ -e "${genotyping_outputs}/${coinf_filename}.fasta" ]; then
                            echo "Coinfection reference database for ${coinf_filename} exists."
                    else
                            echo "making coinf ${coinf_filename} database."
                            grep -A1 "${parent1}" ${lib_prep_database} | grep -v "^--$" > ${genotyping_outputs}/${coinf_filename}.fasta
                            grep -A1 "${parent2}" ${lib_prep_database} | grep -v "^--$" >> ${genotyping_outputs}/${coinf_filename}.fasta
                            [ -n "$parent3" ] && grep -A1 "${parent3}" ${lib_prep_database} | grep -v "^--$" >> ${genotyping_outputs}/${coinf_filename}.fasta
                            [ -n "$parent4" ] && grep -A1 "${parent4}" ${lib_prep_database} | grep -v "^--$" >> ${genotyping_outputs}/${coinf_filename}.fasta

                    fi
                    database_file="${coinf_filename}.fasta"

                    ;;
            *)
                    echo "Unknown type: ${sample_type}"
                    ;;
    esac



    #Usearch files

    # Want to set this up so that it is flexible - right now, there is only one .fastq in every plate/well/ folder,
    # but previously we sorted by segment prior to usearching. we may want to return to doing that
    # shopt -s nullglob means if no files match the string with the wildcard, return nothing, instead of return the wildcard literally
    # theoretically, the glob shouldnt match nothing anyway, but will be helpful if we ever end up


    # Since plate??_well??_${reads_name}.fastq AND a future locus-pre-sorted case (plate??_well??_<locus>_${reads_name}.fastq,
    # with * matching "<locus>_" will run the loop below correctly either 0 times (if there is no fastq for a well),
    # 1 time for the current normal case, or as many times as there are files in the plate/well folder

    shopt -s nullglob
    locus_fastqs=(${cutadapt_outputs}/${plate}/${well}/plate??_well??_*${reads_name}.fastq)
    shopt -u nullglob

    for locus_fastq in "${locus_fastqs[@]}";
    do
            locus_basepath="${locus_fastq##*/}";

            #extract locus name via parameter expansion instead of a sed backreference
            #(the old sed line referenced \1 with no capture group - invalid, silently
            #errored every call, harmless only because $locus was never used downstream)
            locus="${locus_basepath#plate??_well??_}"
            locus="${locus%_${reads_name}.fastq}"
            if [ "$locus" = "${reads_name}.fastq" ]; then
                    locus=""   # no locus token in this filename - current single-file setup (no pre-sorting by segment)
            fi

            echo "reads_name variable = ${reads_name}";
            [ -n "$locus" ] && echo "locus variable = ${locus}";

            #USEARCH STRAIN ASSIGNMENT!!
            #b6 and notmatched fastq filenames will include a locus name only if one is present (like if we add that back in)
            #so since rn the names will fit, this code will prevent loci overwriting each other if we bring back pre-usearch locus sorting
            b6_out="${genotyping_outputs}/${sample_type}/${plate}/${plate}_${well}${locus:+_${locus}}_90_merged.b6"
            notmatched_out="${genotyping_outputs}/${sample_type}/${plate}/${plate}_${well}${locus:+_${locus}}_notmatched.fastq"

            echo
            echo "Running USEARCH on ${locus_basepath} with DB: ${genotyping_outputs}/${database_file} and outputting to: ${b6_out} "
            echo
            usearch -usearch_global "$locus_fastq" -db ${genotyping_outputs}/${database_file} -id 0.90 -blast6out "$b6_out" -strand both -top_hit_only -notmatched "$notmatched_out"

    done

    pars="${database_file%.fasta}"

    fastq="${cutadapt_outputs}/${plate}/${well}/${plate}_${well}_${reads_name}.fastq"
    total_reads=$(($(wc -l < "$fastq") / 4))

    base="${sample_type}_${pars}_${plate}"


    # shopt -s nullglob means if no files match the string with the wildcard, return nothing, instead of return the wildcard literally

    shopt -s nullglob

    #if there are no .b6 files matching the string then the pattern disappears, and sets b6_files to nothing
    #this loop is flexible - sometimes there will be no matching .b6 file, setting nullglob allows that.
    #and in the future if we go bac to pre-sorting by segment we will also ble flexible for that

    #set b6_files to all matching files within the sample_type/plate folders
    b6_files=(${genotyping_outputs}/${sample_type}/${plate}/${plate}_${well}*_90_merged.b6)

    #then unset aka disable the nullglob behavior
    shopt -u nullglob

    #loop thru any .b6 files for the current values of plate_well and use ask to get the total matches for segment and majority strain
    #append them to the summary file:
    for b6_file in "${b6_files[@]}";
    do
            #derive the notmatched reads filename the same way b6_out/notmatched_out were built above
            #I would guess that every single file will also have some reads that dont match,
            #but may not be true for some - for ex, misassigned samples with only a few reads total

            b6_basename="${b6_file##*/}"
            notmatched_file="${b6_file%_90_merged.b6}_notmatched.fastq"

            notmatched_count=0
            #use the standard structure of a fastq file to count the number of reads in the notmatched .fastq file
            if [ -s "$notmatched_file" ]; then
                first_char=$(head -c1 "$notmatched_file")
                if [ "$first_char" = "@" ]; then
                    notmatched_count=$(( $(wc -l < "$notmatched_file") / 4 ))
                elif [ "$first_char" = ">" ]; then
                    notmatched_count=$(grep -c '^>' "$notmatched_file")
                fi
            fi

            awk -v file="${b6_basename}" \
            -v base="${base}" \
            -v outdir="${genotyping_outputs}/strain_assignment_output" \
            -v total_reads="${total_reads}" \
            -v unsorted="${notmatched_count}" '
    {
        # Find the field that looks like CA68_M, PHI6_L, etc. (Strain_Segment) and assign to strain_seg_field
        strain_seg_field = ""
        for (field_index = 1; field_index <= NF; field_index++) {
            if ($field_index ~ /^[A-Za-z0-9]+_[A-Za-z0-9]+$/) {
                strain_seg_field = $field_index
                break
            }
        }

        # If we didnt find a strain_seg_field, hard fail
        if (strain_seg_field == "") {
                printf "ERROR: no strain_segment-like field found in line:\n%s\n", $0 > "/dev/stderr"
            exit 1
        }

        # use the single underscore naming convention (Strain_Segment) to split the name into
        # its two halves: the strain name and the segment name
        split(strain_seg_field, strain_seg_parts, "_")
        segment_name     = strain_seg_parts[2]
        strain_seg_match = strain_seg_field

        reads_for_segment[segment_name]++
        hits_for_strain_seg_match[segment_name, strain_seg_match]++

        # track how many hits each read has (query ID is $1)
        hits_for_read[$1]++
    }

    END {
        summary_file    = outdir "/" base "_strain_assignment_output.txt"
        all_samples_file = outdir "/" base "_strain_assignment_output_all_samples.txt"

        if (length(reads_for_segment) == 0) {
            # no hits at all for this well - summary_file gets nothing (unchanged
            # behavior). all_samples_file gets one explicit row so this well is never
            # silently invisible: segment/match = NONE, all hit counts 0, unsorted
            # = full read count for this well (from -notmatched).
            printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\n",
                   file, "NONE", "NONE", 0, 0, unsorted, 0 >> all_samples_file
        } else {
            # First pass: find the majority (top) match per segment
            # (unsorted is passed in directly from usearch -notmatched output,
            # not derived by subtraction)
            for (segment_name in reads_for_segment) {
                top_match       = ""
                top_match_count = 0

                for (seg_match_key in hits_for_strain_seg_match) {
                    split(seg_match_key, key_parts, SUBSEP)
                    if (key_parts[1] == segment_name && hits_for_strain_seg_match[seg_match_key] > top_match_count) {
                        top_match_count = hits_for_strain_seg_match[seg_match_key]
                        top_match       = key_parts[2]
                    }
                }

                top_count_for_seg[segment_name] = top_match_count
                top_match_for_seg[segment_name] = top_match
            }

            # compute multi-hit reads (reads that hit more than one distinct strain_seg_field)
            multi_hit_reads = 0
            for (read_id in hits_for_read)
                if (hits_for_read[read_id] > 1)
                    multi_hit_reads++

            # second pass: print output for each segment, to both files
            for (segment_name in reads_for_segment) {
                printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\n",
                       file, segment_name, top_match_for_seg[segment_name],
                       top_count_for_seg[segment_name], reads_for_segment[segment_name],
                       unsorted, multi_hit_reads >> summary_file

                printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\n",
                       file, segment_name, top_match_for_seg[segment_name],
                       top_count_for_seg[segment_name], reads_for_segment[segment_name],
                       unsorted, multi_hit_reads >> all_samples_file
            }
        }
    }
        ' "$b6_file"

        echo "DEBUG: base='$base' well='$well' b6_file='$b6_file'"
    done

done

###############################################
# Optional: BLAST the reads usearch couldn't match to anything, against the
# complete library-prep reference database. Runs once at the very end across
# every sample_type/plate/well at once, rather than inside the main loop above,
# so it doesn't reprint the same setup messages hundreds of times.
###############################################

if [[ "$do_blast" == "1" ]]; then

    for tool in blastn makeblastdb seqtk; do
        command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not found on PATH - needed for -b 1 (BLAST notmatched reads). Is your blast conda environment active?"; exit 1; }
    done

    echo "You turned on the option to BLAST the reads usearch couldn't match to the database! Blasting now - this can take a while."

    blast_db_name="${genotyping_outputs}/${lib_prep_database_filename%.*}_blast_nt_db"

    echo "Checking if your complete library prep database file ${lib_prep_database_filename} has already been converted to a BLAST database"
    if [ -e "${blast_db_name}.nin" ] || [ -e "${blast_db_name}.ndb" ]; then
        echo "BLAST database already exists at ${blast_db_name} - skipping makeblastdb"
    else
        echo "Building a new BLAST database from ${lib_prep_database_filename}"
        makeblastdb -in "${lib_prep_database}" -dbtype nucl -out "${blast_db_name}"
    fi

    # matches every notmatched fastq across every sample_type folder at once (coinfection,
    # positive, negative, misassigned), including the optional locus-suffixed filename
    shopt -s nullglob
    notmatched_fastqs=(${genotyping_outputs}/*/plate??/plate??_well??*_notmatched.fastq)
    shopt -u nullglob

    echo "Found ${#notmatched_fastqs[@]} notmatched fastq file(s) to BLAST"

    for notmatched_fastq in "${notmatched_fastqs[@]}";
    do
        [ -s "$notmatched_fastq" ] || continue   # skip empty files - nothing to blast

        fasta="${notmatched_fastq%.fastq}.fasta"
        blast_out="${notmatched_fastq%.fastq}_blastn.tsv"

        seqtk seq -A "$notmatched_fastq" > "$fasta"

        blastn -query "$fasta" -db "${blast_db_name}" -num_threads 60 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid pident length evalue stitle" -out "$blast_out"

        rm -f "$notmatched_fastq"
    done

    echo "FYI - since we had to turn your notmatched fastqs into fastas to BLAST them, we deleted the original fastqs to save some space (the .fasta files and BLAST results are still there)"

elif [[ "$do_blast" == "0" ]]; then
    echo "You turned off BLAST, so we're done!"

else
    echo "You set '\$do_blast' to something besides 0 or 1! I wasn't sure what to do, so I skipped BLASTing this time :("
fi

echo
echo "completed genotyping with usearch - look for your results in genotyping_outputs/strain_assignment_output"
