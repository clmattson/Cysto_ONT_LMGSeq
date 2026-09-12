#!/bin/bash

#current run:
#bash demultiplexing.sh -d /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall -r /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall/FBC73506_fastq_pass_9ee158db_6bd36f99_0.fastq -p /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall/plate_barcodes.fasta -w /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall/well_barcodes.fasta -c 0 -l 100

# ############## Section I: collect user inputs, move existing data #################################

#set up function to get creation timestamp for existing directories so that if this script is run, old data wont get overwritten
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

# function to print a message followed by a blank line for readability
printf_nl() {
    printf '%s\n\n' "$*"
}

# quick fail-fast check so a forgotten conda activation shows up immediately with a
# clear message, instead of a confusing error several steps into the run
for tool in porechop cutadapt seqkit; do
        command -v "$tool" >/dev/null || { echo "Oops! '$tool' not found on PATH - did you activate the right conda environment (cutadapt)?"; exit 1; }
done

#main script, input to gather:
echo "flag info:"
echo " -d = path to the desired working dir"
echo " -r = path to reads"
echo " -p = <plate_barcode.fasta> file with plate barcodes fasta with path"
echo " -w = <well_barcode.fasta> file with well barcodes fasta with path"
echo " -c = do_porechop: 0 or 1 for whether or not to re-do the porechop step"
echo " -l = min_length: minimum desired length cutoff for filtering reads"

#do_porechop is 0 for no porechop, like if its already done 1 for yes porechop. this is for debugging by courtney

#initialize user-input objects
demuxed_path=''
reads=''
plate_barcodes=''
well_barcodes=''
min_length=''

print_usage() {
  printf "Usage: ..."
}

#read in user-input flags, defined above
while getopts d:r:p:w:c:l: flag
do
    case "${flag}" in
        d) working_dir=${OPTARG};;
        r) reads_path=${OPTARG};;
        p) plate_barcodes=${OPTARG};;
        w) well_barcodes=${OPTARG};;
        c) do_porechop=${OPTARG};;
        l) min_length=${OPTARG};;
    esac
done

#dont forget to load and activate your conda env!
#module load conda; conda activate cutadapt
# Note: Start analysis with basecalled data from Dorado - run dorado with trimming DISABLED
# Reads should be in one large .fastq file

echo
printf "Your current flag inputs:\n -d = ${working_dir}\n -r = ${reads_path}\n p = ${plate_barcodes}\n w = ${well_barcodes}\n -c = ${do_porechop}\n -l = ${min_length}"
echo

# MOVE EXISTING OUTPUTS and RE-SET

# move old outputs to timestamped folder so they won't be overwritten
if [ -d "${working_dir}/cutadapt_outputs" ]; then
    old_dir_cut="${working_dir}/cutadapt_outputs"

    # get earliest timestamp (mtime/ctime)
    ts_cut=$(get_dir_timestamp "$old_dir_cut")

    echo "a cutadapt_outputs directory already exists in the location you set - moving existing cutadapt_outputs directory to cutadapt_outputs_from_${ts_cut} "
        echo
    mv "$old_dir_cut" "${old_dir_cut}_from_${ts_cut}"
fi

reads_name="${reads_path%.*}";
reads_name="${reads_name##*/}";
printf_nl "reads_name variable = ${reads_name}"
printf_nl "We will now use a combination of porechop and cutadapt to demulitplex the data:"

# save the real terminal file descriptors so each section below can redirect its own
# output into its own log file, then restore back to the terminal afterward - this
# keeps porechop's log and cutadapt's log separate, and means the person running this
# doesn't have to remember to add their own tee command to log when they run the script
exec 3>&1 4>&2

# ################## Section II: Porechop ############################################3

if [[ "$do_porechop" == "1" ]]; then

        # move old porechop outputs to time stamped dir - this check has to happen
        # BEFORE mkdir -p creates a fresh porechop_outputs below, otherwise the
        # directory would always already exist by the time we check for it
        if [ -d "${working_dir}/porechop_outputs" ]; then
            old_dir_pore="${working_dir}/porechop_outputs"
            ts_pore=$(get_dir_timestamp "$old_dir_pore")
            printf_nl "FYI, a porechop_outputs dir already exists - moving existing porechop_outputs directory to porechop_outputs_from_${ts_pore}"
            mv "$old_dir_pore" "${old_dir_pore}_from_${ts_pore}"
        fi

        #make directory for new porechop outputs
        mkdir -p ${working_dir}/porechop_outputs
        porechop_outputs="${working_dir}/porechop_outputs"

        porechop_log="${porechop_outputs}/porechop_screen.log"
        exec > >(tee -a "$porechop_log") 2>&1

        printf_nl "You turned on porechop! Full porechop-stage output is being saved to: ${porechop_log}"
        printf_nl "made new directory for porechop outputs"

        #Note: the next section of code breaks the input fastq's into chunks for porechopping
        # This is necessary because Porechop loads the entire input file into memory, which takes an enormous amt of memory for a large fastq file. It is not fixed by --threads
        # This is a known porechop bug, but unfortunately the tool is no longer maintained: https://github.com/rrwick/Porechop/issues/77
        # It worked fine for me with 5 sections, but the number could be increased if needed.
        # We also use declare and wait to do some processing simultaneously on porechopped reads (calculate some stats before and after porechopping and re-naming duplicate reads names on chopped reads)
        # FYI - duplicate read names is a known and corrected bug in porechop_ABI, but since porechop original is not maintained it is not fixed there: https://github.com/nf-core/mag/issues/840

        # Number of parts to split the fastq into - 5 for now. Get number of reads to add to each chunk
        num_fastq_sections=5
        fastq_total_reads=$(wc -l ${reads_path} | awk '{print $1 / 4}')
        lines_per_section=$(( (fastq_total_reads / num_fastq_sections) * 4 ))

        split -l "${lines_per_section}" -d --additional-suffix=.fastq \
            "${reads_path}" "${working_dir}/${reads_name}_chunk"

        # track background post-processing jobs
        declare -a post_jobs=()

        #run porechop on each fastq chunk
        for chunk in 00 01 02 03 04 05; do

            infile="${working_dir}/${reads_name}_chunk${chunk}.fastq"
            outfile="${porechop_outputs}/${reads_name}_porechop_chunk${chunk}.fastq"
            logfile="${porechop_outputs}/${reads_name}_porechop_chunk${chunk}.log"
            fixed_out="${porechop_outputs}/${reads_name}_porechop_chunk${chunk}_unique.fastq"
            stats_file="${porechop_outputs}/${reads_name}_porechop_chunk${chunk}_stats.txt"

            printf_nl "Running porechop on chunk ${chunk}"

            # STEP 1: PORECHOP RUNS - change any porechop settings here!
            porechop -i "$infile" --verbosity 2 --end_threshold 70 --middle_threshold 80 \
                --extra_end_trim 0 --end_size 150 --min_split_read_size 200 \
                --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 \
                --min_trim_size 8 -o "$outfile" > "$logfile"

            printf_nl "Finished porechop chunk ${chunk}"

            # POST-PROCESSING RUNS IN BACKGROUND - speed things up and get info on poreshopped data
            (

                printf_nl "===== Porechop info for chopping of chunk ${chunk}; from log file: $logfile ====="
                head -n 22 "$logfile"
                echo
                tail -n +23 "$logfile" | grep -F "adapters"
                echo

                # read stats BEFORE - calc total reads and % of reads longer than 1000 bp
                #num reads
                total_before=$(($(wc -l < "$infile") / 4))
                #reads over 1000 bp
                long_before=$(awk 'NR%4==2 { if(length($0) > 1000) c++ } END { print c+0 }' "$infile")
                #percent of total reads that were longer than 1000 bp
                pct_before=$(awk -v a="$long_before" -v b="$total_before" 'BEGIN { printf("%.2f", (a/b)*100) }')

                # read stats AFTER. same as above - calculate % of reads longer than 1000% after porechopping
                total_after=$(grep -c "^@" "$outfile")
                long_after=$(awk 'NR%4==2 { if(length($0) > 1000) c++ } END { print c+0 }' "$outfile")
                pct_after=$(awk -v a="$long_after" -v b="$total_after" 'BEGIN { printf("%.2f", (a/b)*100) }')

                echo "++++++ Read-length summary for chunk ${chunk} (what % is >1 kb?) +++++++"
                echo "  Before Porechop: ${pct_before}% = (${long_before} / ${total_before}) "
                printf_nl "  After  Porechop: ${pct_after}% = (${long_after} / ${total_after}) "

                # save this chunk's raw numbers so the parent shell can sum them into a
                # combined summary later - background subshells can't set variables that
                # persist back in the parent shell, so a small file is the simplest handoff
                {
                    echo "total_before=${total_before}"
                    echo "total_after=${total_after}"
                    echo "long_before=${long_before}"
                    echo "long_after=${long_after}"
                } > "$stats_file"

                echo "Porechopping can leave reads behind with identical names."
                printf_nl "Now fixing duplicate read names for chunk ${chunk} "
                #use awk to appeand readname_<num> to read names to eliminate identical read names after porechopping
                awk '
                BEGIN { FS=" "; OFS=" " }
                NR%4==1 {
                    id=$1
                    count[id]++
                    if (count[id] > 1) {
                        $1 = id "_" count[id]
                    }
                }
                { print }
                ' "$outfile" > "$fixed_out"
                echo
                echo "Finished fixing names for chunk ${chunk} "
            ) &

            # record background job PID
            post_jobs+=($!)

        done
	echo
	printf_nl "All porechop chunks complete!"
        printf_nl "Waiting for post-processing jobs to finish..."
        # wait for all background post-processing jobs
        for pid in "${post_jobs[@]}"; do
            wait "$pid"
        done

        ###############################################
        # Combined porechop summary across all chunks
        ###############################################
        printf_nl "===== COMBINED porechop read-length summary across all ${num_fastq_sections} chunks ====="

        combined_total_before=0
        combined_total_after=0
        combined_long_before=0
        combined_long_after=0

        combined_start_trim_reads=0
        combined_start_trim_bp=0
        combined_end_trim_reads=0
        combined_end_trim_bp=0
        combined_middle_split_reads=0

        for stats_file in ${porechop_outputs}/${reads_name}_porechop_chunk*_stats.txt; do
            # each stats file just contains var=value lines we wrote ourselves above,
            # so sourcing it directly is safe here
            source "$stats_file"
            combined_total_before=$(( combined_total_before + total_before ))
            combined_total_after=$(( combined_total_after + total_after ))
            combined_long_before=$(( combined_long_before + long_before ))
            combined_long_after=$(( combined_long_after + long_after ))

            # the matching chunk .log file (porechop's own raw tool output) is what has
            # the adapter-trim lines - derive its path from the stats filename rather
            # than tracking a separate variable
            chunk_log="${stats_file%_stats.txt}.log"

            # each of these lines may or may not be present (e.g. a chunk where porechop
            # found no adapters at all, like a 1-read chunk, prints "No adapters found"
            # instead) - default cleanly to 0 rather than erroring if grep finds nothing
            line_start=$(grep "reads had adapters trimmed from their start" "$chunk_log")
            if [ -n "$line_start" ]; then
                start_reads=$(echo "$line_start" | awk -F' / ' '{print $1}' | tr -d ',')
                start_bp=$(echo "$line_start" | grep -oE '[0-9,]+ bp removed' | grep -oE '[0-9,]+' | tr -d ',')
            else
                start_reads=0
                start_bp=0
            fi
            start_reads=${start_reads:-0}
            start_bp=${start_bp:-0}

            line_end=$(grep "reads had adapters trimmed from their end" "$chunk_log")
            if [ -n "$line_end" ]; then
                end_reads=$(echo "$line_end" | awk -F' / ' '{print $1}' | tr -d ',')
                end_bp=$(echo "$line_end" | grep -oE '[0-9,]+ bp removed' | grep -oE '[0-9,]+' | tr -d ',')
            else
                end_reads=0
                end_bp=0
            fi
            end_reads=${end_reads:-0}
            end_bp=${end_bp:-0}

            line_middle=$(grep "reads were split based on middle adapters" "$chunk_log")
            if [ -n "$line_middle" ]; then
                middle_reads=$(echo "$line_middle" | awk -F' / ' '{print $1}' | tr -d ',')
            else
                middle_reads=0
            fi
            middle_reads=${middle_reads:-0}

            combined_start_trim_reads=$(( combined_start_trim_reads + start_reads ))
            combined_start_trim_bp=$(( combined_start_trim_bp + start_bp ))
            combined_end_trim_reads=$(( combined_end_trim_reads + end_reads ))
            combined_end_trim_bp=$(( combined_end_trim_bp + end_bp ))
            combined_middle_split_reads=$(( combined_middle_split_reads + middle_reads ))
        done

        combined_pct_before=$(awk -v a="$combined_long_before" -v b="$combined_total_before" 'BEGIN { printf("%.2f", (a/b)*100) }')
        combined_pct_after=$(awk -v a="$combined_long_after" -v b="$combined_total_after" 'BEGIN { printf("%.2f", (a/b)*100) }')

	echo
        echo "  TOTAL reads before porechop: ${combined_total_before}"
        echo "  TOTAL reads after porechop:  ${combined_total_after}"
        echo "  Before Porechop: ${combined_pct_before}% >1kb (${combined_long_before} / ${combined_total_before})"
        echo "  After  Porechop: ${combined_pct_after}% >1kb (${combined_long_after} / ${combined_total_after})"
        echo
        echo "  TOTAL reads with adapters trimmed from start: ${combined_start_trim_reads} (${combined_start_trim_bp} bp removed)"
        echo "  TOTAL reads with adapters trimmed from end:   ${combined_end_trim_reads} (${combined_end_trim_bp} bp removed)"
        printf_nl "  TOTAL reads split on middle adapters:         ${combined_middle_split_reads}"

        printf_nl "Now pasting all porechopped chunks back into one file: ${porechop_outputs}/${reads_name}_porechop.fastq"
        cat ${porechop_outputs}/${reads_name}_porechop_chunk*_unique.fastq \
            > ${porechop_outputs}/${reads_name}_porechop.fastq

        ###############################################
        # Clean up chunk-level intermediate files now that the merged file exists.
        # Everything deleted here is fully redundant with ${reads_name}_porechop.fastq -
        # nothing downstream reads these again. Kept: the per-chunk .log files (porechop's
        # own tool output, small, useful reference) and the per-chunk _stats.txt files
        # (tiny, already summarized into the combined summary above but nice as a raw record).
        ###############################################
        printf_nl "Cleaning up chunk-level intermediate fastqs (raw split chunks, pre-rename porechop output, and per-chunk renamed output) - all fully redundant now that ${reads_name}_porechop.fastq exists"
        rm -f ${working_dir}/${reads_name}_chunk??.fastq
        rm -f ${porechop_outputs}/${reads_name}_porechop_chunk??.fastq
        rm -f ${porechop_outputs}/${reads_name}_porechop_chunk??_unique.fastq

        printf_nl "Porechop complete: executed porechop (original) for splitting reads on landing pads"

        # restore stdout/stderr to the terminal before leaving this section
        exec 1>&3 2>&4

elif [[ "$do_porechop" == "0" ]]; then
        echo "You turned off porechop, so we will use the existing, already-chopped reads in the porechop_outputs directory"

        if [ ! -d "${working_dir}/porechop_outputs" ]; then
        echo "silly, you turned off porechop but the ${working_dir}/porechop_outputs doesnt exist. you cant do that!"
        exit 1
        fi

        porechop_outputs="${working_dir}/porechop_outputs"
else
        echo "You set '$do_porechop' to something besides 0 or 1, exiting the enitre script! Pick one and run again." >&2
        exit 1
fi

# ########################## Section III: Demultiplexing w Cutadapt ####################################

#make directory for cutadapt outputs:
mkdir -p ${working_dir}/cutadapt_outputs
cutadapt_outputs="${working_dir}/cutadapt_outputs"

cutadapt_log="${cutadapt_outputs}/cutadapt_screen.log"
exec > >(tee -a "$cutadapt_log") 2>&1
printf_nl "Full cutadapt-stage script output is being saved to: ${cutadapt_log}"
printf_nl "(cutadapt itself also writes its own per-plate/per-well .log files inside cutadapt_outputs - this file is the script's own narration, not a replacement for those)"

echo "Is the cutadapt_outputs directory variable set correctly? ~~~~cutadapt_outputs='$cutadapt_outputs'~~ and does it contain 0 things as it should? :"
ls -l "$cutadapt_outputs" | wc -l
echo

#DEMULTIPLEXING - STEP 2 - PLATE
echo "directory made and check complete - now doing plate-level demultiplexing"
#use plate_barcodes.fasta file to search and demultiplex PLATE barcodes with cutadapt. Higher -O
#07-30-2026 CM update: add length filtering - user-specified length value
# --times 10: re-search each read for more adapters after trimming one, so info-file shows if a read had >1 adapter; -j 0: use all cores
cutadapt -a file:${plate_barcodes} -O 14 --revcomp -e 0.15 --times 10 --cores=0 --minimum-length ${min_length} --too-short-output ${cutadapt_outputs}/${reads_name}_cutadapt_under_${min_length}_bp.fastq --info-file ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop_INFO.tsv -o ${cutadapt_outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${porechop_outputs}/${reads_name}_porechop.fastq > ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop.log

printf_nl "completed demultiplexing step 2 - cutadapt plate identification!"

# Use the INFO .tsv file from plate demuxing with cutadapt to search for any reads that matched to more than one plate, and filter them out
#plate_info="${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop_INFO.tsv"
plate_multi_hit_pairs="${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop_multi_adapter_pairs.txt"
plate_multi_hit_ids="${cutadapt_outputs}/${reads_name}_multi_adapter_ids.txt"

#at least right now, I only want to filter out read ID's that match to multiple DIFFERENT barcodes with cutadapt,
#but leave those that match >1x to the same barcode. So here we fetch the nonunique read names from the info file and check if they have the same match
#easiest way will bepull out duplicated read names in the tsv, and then search for the barcode names that we just used to run cutadapt to see if they are unique or different

plate_names=$(grep '^>' "$plate_barcodes" | sed 's/^>//')

awk -v names="$plate_names" '
BEGIN {
    #add the plate names to an array, separated by newline chars
    split(names, plate_names_array, "\n")
    #create a lookup table of barcode names
    for (i in plate_names_array) valid_plates[plate_names_array[i]] = 1
}
FNR==NR { if (FNR%4==1) { sub(/^@/,"",$1); too_short[$1]=1 }; next }
$1 in too_short { next }
{
   #awk automatically loops thru all lines in the file
   #within each line, loop thru every field up to the tot items in the line, uses awk $2 = $2nd item syntax
    for (field = 2; field <= NF; field++)
        #if the item in the line is in the plate barcode array, then print the first item in the line as well as the matched column
        if ($field in valid_plates) {
            print $1, $field
            break
        }
}
' "${cutadapt_outputs}/${reads_name}_cutadapt_under_${min_length}_bp.fastq" "${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop_INFO.tsv" | sort -u > "$plate_multi_hit_pairs"
#get just the list of unique read/barcode hit PAIRS (1st col only can contain dups)

#save non-unique values from the first field of that list (the read ids), so reads that by def would have matched more than once, but to diff barcodes
awk '{print $1}' "$plate_multi_hit_pairs" | uniq -d > "$plate_multi_hit_ids"

#get the list of duplicated reads and the barcodes they matched
awk 'NR==FNR { flagged[$1]; next } $1 in flagged' "$plate_multi_hit_ids" "$plate_multi_hit_pairs" > "${plate_multi_hit_pairs}.tmp" \
  && mv "${plate_multi_hit_pairs}.tmp" "$plate_multi_hit_pairs"

num_plate_multi=$(wc -l < "$plate_multi_hit_ids")
printf_nl "We found ${num_plate_multi} reads with >1 distinct plate adapter (excluding those already removed for being too short). Now lets remove them"

for f in ${cutadapt_outputs}/plate??_${reads_name}_cutadapt_porechop.fastq; do
     # strip out flagged multi-adapter reads before the empty check below
    if [ -s "$plate_multi_hit_ids" ]; then  #check if there are 1 or more reads flagged as hitting to multiple adapters
        multi_plate_filtered="${f%.fastq}_filtered.fastq"
        #use seqkit to reverse grep (aka filter) any of the multi hit reads ids
        seqkit grep -v -f "$plate_multi_hit_ids" "$f" -o "$multi_plate_filtered"
        #repplace the original fastq qith the filtered one
        mv "$multi_plate_filtered" "$f"
        echo "filtered reads with multiple plate hits for ${f}"
    else
        echo "no reads were found to have matched to multiple plate barcodes for ${f}"
    fi

    #delete any empty plate barcodes.
    lines=$(wc -l < "$f")
    if [ "$lines" -lt 4 ]; then
        rm -f "$f"
    fi
done

#Use find to loop thru the plate-demuxed files we just created
find "${cutadapt_outputs}" -type f -name 'plate??_*.fastq' | while read -r plate_file_path; do
        printf_nl "entered well cutadapt loop, plate_file_path = ${plate_file_path}"

        #Move the plate_ demultiplexed files we just made into directories based off the file names:
        plate_file_name="$(basename "$plate_file_path")"
        plate="${plate_file_name%%_*}"
        plate_dir="${cutadapt_outputs}/${plate}";
        mkdir -p "${plate_dir}"
        mv "${plate_file_path}" "${plate_dir}/";
        printf_nl "File sorted into directory 'plate_dir' : ${plate_dir}"

        #Ok we want to execute the well-demultiplexing step once for each plate file. so include it in this loop:
        echo "we just demuxed cutadapt_outputs by PLATE, and made a new cutadapt_outputs/plate?? directory for each plate barcode. We are about to start demultiplexing by well with cutadapt.\n"
        echo "lets check that ${plate_dir} which we just set, exist and contains the correct number of items with ls"
        ls -1 "$plate_dir" | wc -l

        #DEMULTIPLEXING - STEP 3 - WELL
        #07-30-2026 CM update: change to O=18 (75%)
        #07-30-2026 CM update: change the 5' adapter search to cut from the RIGHTMOST match if multiple matches to the SAME 5' adapter are found in a single read
        printf_nl "executing demultiplexing step 3: cutadapt search for well barcodes! input file=${plate_dir}/${plate_file_name}"

        #70-31-2026 CM update: use an array to expand all the barcode names in well barcodes file because file: and ;rightmost arent compatible in cutadapt
        # NEW: work around a cutadapt bug (reproduced on 5.2 and current dev build) where
        # file:...;rightmost together crash with "unexpected keyword argument 'rightmost'".
        # Build individual -g "name=SEQ;rightmost" args from the fasta instead of using file:.
        declare -a well_g_args=()
        while read -r wb_name && read -r wb_seq; do
                wb_name="${wb_name#>}"
                well_g_args+=(-g "${wb_name}=${wb_seq};rightmost")
        done < "${well_barcodes}"

        cutadapt "${well_g_args[@]}" -O 18 --revcomp -e 0.15 --times 10 --cores=0 --info-file ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_dir}/${plate_file_name} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log

    #Removing reads that were matched to multiple wells:
        # Use the INFO .tsv file from well demuxing with cutadapt to search for any reads that matched to more than one well, and filter them out
        well_info="${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop_INFO.tsv"
        well_multi_hit_pairs="${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop_multi_adapter_pairs.txt"
        well_multi_hit_ids="${plate_dir}/${plate}_well_multi_adapter_ids.txt"

        #same logic as the plate-level filter above: only filter read IDs that match multiple DIFFERENT well barcodes,
        #leave reads that match >1x to the same barcode alone

        well_names=$(grep '^>' "$well_barcodes" | sed 's/^>//')

        awk -v names="$well_names" '
        BEGIN {
            #add the well names to an array, separated by newline chars
            split(names, well_names_array, "\n")
            #create a lookup table of barcode names
            for (i in well_names_array) valid_wells[well_names_array[i]] = 1
        }
        {
           #awk automatically loops thru all lines in the file
           #within each line, loop thru every field up to the tot items in the line, uses awk $2 = $2nd item syntax
            for (field = 2; field <= NF; field++)
                #if the item in the line is in the well barcode array, then print the first item in the line as well as the matched column
                if ($field in valid_wells) {
                    print $1, $field
                    break
                }
        }
        ' "$well_info" | sort -u > "$well_multi_hit_pairs" #get full list of pairs, minus any FULL duplicates (one read matches to same barcode >1x
        awk '{print $1}' "$well_multi_hit_pairs" | uniq -d > "$well_multi_hit_ids" #save non-unique values from the first field of that list (the read ids), so reads that by def would have matched more than once, but to diff barcodes
        #get the list of duplicated reads and the barcodes they matched
        awk 'NR==FNR { flagged[$1]; next } $1 in flagged' "$well_multi_hit_ids" "$well_multi_hit_pairs" > "${well_multi_hit_pairs}.tmp" \
  && mv "${well_multi_hit_pairs}.tmp" "$well_multi_hit_pairs"

        num_well_multi=$(wc -l < "$well_multi_hit_ids")
        printf_nl "We found ${num_well_multi} reads with >1 distinct well adapter under (${plate}). Now lets remove them"

        for f in "${plate_dir}"/"${plate}"_well??_${reads_name}_cutadapt_porechop.fastq; do
             # strip out flagged multi-adapter reads before the empty check below
            if [ -s "$well_multi_hit_ids" ]; then  #check if there are 1 or more reads flagged as hitting to multiple adapters
                multi_well_filtered="${f%.fastq}_filtered.fastq"
                #use seqkit to reverse grep (aka filter) any of the multi hit read ids
                seqkit grep -v -f "$well_multi_hit_ids" "$f" -o "$multi_well_filtered"
                #replace the original fastq with the filtered one
                mv "$multi_well_filtered" "$f"
                echo "filtered reads with multiple well hits for ${f}"
            else
                echo "no reads were found to have matched to multiple well barcodes for ${f}"
            fi

            #delete any empty well fastqs
            lines=$(wc -l < "$f")
            if [ "$lines" -lt 4 ]; then
                rm -f "$f"
            fi
        done

        #move plate??_well??_ demultiplexed files to well folders:
        #plate_dir is updated for each value of plate
        # Make sure it's a directory
        [ -d "$plate_dir" ] || continue

        #Loop through matching files inside the plate directory
        for plate_well_file_path in "${plate_dir}"/"${plate}"_well??_*.fastq; do
                # Extract the filename
                plate_well_file_name="$(basename "$plate_well_file_path")";

                # Extract the well ID (e.g., well01)
                well="${plate_well_file_name#*_}";
                well="${well%%_*}";

                # Create the well subdirectory
                plate_well_dir="$plate_dir/$well";
                mkdir -p "$plate_well_dir";

                # Move the file into the well subdirectory
                mv "${plate_well_file_path}" "${plate_well_dir}/";

                #while still looping thru values of plate and well, do cutadapt search for primers

                #DEMULTIPLEXING - STEP 4 - SEGMENT - TURN OFF FOR NOW
                #okay demultiplex by plaque, input = plate-demuxed files; -O is smaller bc the primers are shorter
                #cutadapt -a small=CTTTCGTACAACCGAGTAGG...CTCCTGAAGTATCTCACGCC -a medium=CGCTACGGCGGTATTGTC...GCTCACCAAGTAAGGTGTAGTAT -a large=TCGATGTTCAACTACTACGC...GCGAGACTCGCTTTGC -O 10 --revcomp -e 0.15 --cores=0 --info-file ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_well_dir}/${plate}_${well}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_well_dir}/${plate_well_file_name} > ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop.log;

         # If directory is now empty, delete it
         rmdir "${plate_well_dir}" 2>/dev/null
       #close final re-organization loop
       done

done

# restore stdout/stderr to the terminal before the final summary message
exec 1>&3 2>&4

printf_nl "Dont forget!! I moved your previous cutadapt_outputs and porechop_outputs directories to cutadapt_outputs_from_${ts_cut} & porechop_outputs_from_${ts_pore} :)"
echo "Porechop-stage log saved to: ${porechop_outputs}/porechop_screen.log"
echo "Cutadapt-stage log saved to: ${cutadapt_outputs}/cutadapt_screen.log"
