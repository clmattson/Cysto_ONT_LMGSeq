#!/bin/bash


# NEW: without this, an unmatched glob (e.g. no files matching plate??_well??_*.fastq)
# doesn't just skip - bash leaves the literal pattern string (with the ?? still in it)
# and loops/commands downstream try to operate on that nonexistent literal filename.
# nullglob makes non-matching globs expand to nothing instead.
#shopt -s nullglob

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

echo "flag info:"
echo " -d = path to the desired working dir"
echo " -r = path to reads"
echo " -p = <plate_barcode.fasta> file with plate barcodes fasta with path"
echo " -w = <well_barcode.fasta> file with well barcodes fasta with path"
echo " -c = do_porechop: 0 or 1 for whether or not to re-do the porechop step"
echo " -l = min_length: minimum desired length cutoff for filtering reads"

#do_porechop is 0 for no porechop, like if its already done 1 for yes porechop. this is for debugging by courtney

demuxed_path=''
reads=''
plate_barcodes=''
well_barcodes=''
min_length=''



print_usage() {
  printf "Usage: ..."
}

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


#module load conda
#conda activate cutadapt

echo
printf "Your current flag inputs:\n -d = ${working_dir}\n -r = ${reads_path}\n p = ${plate_barcodes}\n w = ${well_barcodes}\n -c = ${do_porechop}\n -l = ${min_length}"
echo

# MOVE EXISTING OUTPUTS and RE-SET

# move old outputs to timestamped folder so they won't be overwritten
if [ -d "${working_dir}/cutadapt_outputs" ]; then
    old_dir_cut="${working_dir}/cutadapt_outputs"

    # get earliest timestamp (mtime/ctime)
    ts_cut=$(get_dir_timestamp "$old_dir_cut")

    echo "moving cutadapt_outputs directory to cutadapt_outputs_${ts_cut} "
        echo
    mv "$old_dir_cut" "${old_dir_cut}_from_${ts_cut}"
fi

# Start analysis with basecalled data from Dorado - run dorado with trimming DISABLED
# Reads should be in one large Fastq


reads_name="${reads_path%.*}";
reads_name="${reads_name##*/}";
echo "reads_name variable = ${reads_name} ";
echo
echo "We will now use a combination of porechop and cutadapt to demulitplex the data:"
echo

if [[ "$do_porechop" == "1" ]]; then

          # move old porechop outputs to time stamped dir
        if [ -d "${working_dir}/porechop_outputs" ]; then
            old_dir_pore="${working_dir}/porechop_outputs"
            ts_pore=$(get_dir_timestamp "$old_dir_pore")
            echo "moving porechop_outputs directory to porechop_outputs_${ts_pore}"
            echo
        mv "$old_dir_pore" "${old_dir_pore}_from_${ts_pore}"
        fi


        #make directory for new porechop outputs:
        mkdir -p ${working_dir}/porechop_outputs
        porechop_outputs="${working_dir}/porechop_outputs"

        echo "made new directoy for porechop outputs"
        echo

 #number of parts to split the fastq into
        num_fastq_sections=5
        fastq_total_reads=$(wc -l ${reads_path} | awk '{print $1 / 4}')
        lines_per_section=$(( (fastq_total_reads / num_fastq_sections) * 4 ))

        split -l "${lines_per_section}" -d --additional-suffix=.fastq \
            "${reads_path}" "${working_dir}/${reads_name}_chunk"

        # track background post-processing jobs

        declare -a post_jobs=()

        for chunk in 00 01 02 03 04 05; do

            infile="${working_dir}/${reads_name}_chunk${chunk}.fastq"
            outfile="${porechop_outputs}/${reads_name}_porechop_chunk${chunk}.fastq"
            logfile="${porechop_outputs}/${reads_name}_porechop_chunk${chunk}.log"
            fixed_out="${porechop_outputs}/${reads_name}_porechop_chunk${chunk}_unique.fastq"

            echo
            echo "Running porechop on chunk ${chunk}"

                echo
            # PORECHOP RUNS SERIAL (SAFE)
            porechop -i "$infile" --verbosity 2 --end_threshold 70 --middle_threshold 80 \
                --extra_end_trim 0 --end_size 150 --min_split_read_size 200 \
                --extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0 \
                --min_trim_size 8 -o "$outfile" > "$logfile"

            echo "Finished porechop chunk ${chunk}"

            ###############################################
            # POST-PROCESSING RUNS IN BACKGROUND (FAST)
            ###############################################
            (
                echo
                echo "===== Porechop info for chopping of chunk ${chunk}; from log file: $logfile ====="
                echo

                head -n 22 "$logfile"
                echo

                tail -n +23 "$logfile" | grep -F "adapters"
                echo
                echo

                # read stats BEFORE
                total_before=$(($(wc -l < "$infile") / 4))
                long_before=$(awk 'NR%4==2 { if(length($0) > 1000) c++ } END { print c+0 }' "$infile")
                pct_before=$(awk -v a="$long_before" -v b="$total_before" 'BEGIN { printf("%.2f", (a/b)*100) }')

                # read stats AFTER
                total_after=$(grep -c "^@" "$outfile")
                long_after=$(awk 'NR%4==2 { if(length($0) > 1000) c++ } END { print c+0 }' "$outfile")
                pct_after=$(awk -v a="$long_after" -v b="$total_after" 'BEGIN { printf("%.2f", (a/b)*100) }')

                echo "++++++ Read-length summary for chunk ${chunk} (what % is >1 kb?) +++++++"
                echo "  Before Porechop: ${pct_before}%  (${long_before} / ${total_before}) "
                echo "  After  Porechop: ${pct_after}%  (${long_after} / ${total_after}) "

                echo
                echo "Now fixing duplicate read names for chunk ${chunk} "

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
        echo "All porechop chunks complete!"
        echo "Waiting for post-processing jobs to finish..."

        echo
        # wait for all background post-processing jobs
        for pid in "${post_jobs[@]}"; do
            wait "$pid"
        done

        echo "Now pasting all porechopped chunks back into one file: ${porechop_outputs}/${reads_name}_porechop.fastq"

        echo
        cat ${porechop_outputs}/${reads_name}_porechop_chunk*_unique.fastq \
            > ${porechop_outputs}/${reads_name}_porechop.fastq

        echo "executed porechop for splitting reads on landing pads"
        echo

elif [[ "$do_porechop" == "0" ]]; then
        echo "You turned off porechop so we will use the existing, already-chopped reads in the porechop_outputs directory"

        if [ ! -d "${working_dir}/porechop_outputs" ]; then
        echo "silly, you turned off porechop but the ${working_dir}/porechop_outputs doesnt exist. you cant do that!"
        exit 1
        fi

        porechop_outputs="${working_dir}/porechop_outputs"
else
        echo "You set '$do_porechop' to something besides 0 or 1, now i crashhhh :( " >&2
        exit 1
fi


#Ok lets make directories for the demuliplexing output:
#this makes a folder for every plate barcode and a subfolder for every well barcode. could later change to do based on the files output by cutadapt instead
#for plate_dir in $(grep '^>' ${plate_barcodes} | sed 's/^>//');
#do mkdir "${plate_dir}";
#       for well_dir in $(grep '^>' ${well_barcodes} | sed 's/^>//');
#       do mkdir "${plate_dir}/${well_dir}";
        #done
#done

#make directory for cutadapt outputs:
mkdir -p ${working_dir}/cutadapt_outputs
cutadapt_outputs="${working_dir}/cutadapt_outputs"

echo "Is the cutadapt_outputs directory variable set correctly? ~~~~cutadapt_outputs='$cutadapt_outputs'~~ and does it contain 0 things as it should? :"
echo

ls -l "$cutadapt_outputs" | wc -l

#DEMULTIPLEXING - STEP 2 - PLATE
#use plate_barcodes.fasta file to search and demultiplex PLATE barcodes with cutadapt. Higher -O
#07-30-2026 CM update: add length filtering - user-specified length value
# --times 10: re-search each read for more adapters after trimming one, so info-file shows if a read had >1 adapter; -j 0: use all cores
cutadapt -a file:${plate_barcodes} -O 14 --revcomp -e 0.15 --times 10 --cores=0 --minimum-length ${min_length} --info-file ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop_INFO.tsv -o ${cutadapt_outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${porechop_outputs}/${reads_name}_porechop.fastq > ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop.log
#cutadapt -a file:${plate_barcodes} -O 14 --action=lowercase --revcomp -e 0.15 -o ${cutadapt_outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${working_dir}/${reads_name}_porechop.fastq > ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop.log

echo "completed demultiplexing step 2 - cutadapt plate identification!"
echo

########################################################################
# NEW: ADAPTER QC - keep reads with exactly 1 distinct plate adapter
# (any repeat count); discard reads with 0 or >1 distinct adapters.
# Uses the info-file the command above already made, no extra cutadapt run.
########################################################################

qc_info="${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop_INFO.tsv"
multi_ids="${cutadapt_outputs}/${reads_name}_multi_adapter_ids.txt"

awk -F'\t' '$(NF-1) != "-1" {          # skip no-match rows (they end in ...  -1  <seq>)
    n = split($1, a, /[ \t]/)          # read id = first whitespace token (ONT header tags ignored)
    print a[1], $(NF-3)                # adapter name is always 3 fields from the end
}' "$qc_info" \
  | sort -u \
  | awk '{print $1}' \
  | uniq -d \
  > "$multi_ids"
# sort -u dedupes (id,adapter) pairs, collapsing same-adapter repeats
# uniq -d: id still appears >1x after that dedupe = >1 DISTINCT adapter
# no-match reads need no handling here - they never land in a plate??_*.fastq
# file to begin with (cutadapt routes them to its own "unknown" output)

n_multi=$(wc -l < "$multi_ids" 2>/dev/null || echo 0)
echo "Adapter QC: flagged ${n_multi} reads with >1 distinct plate adapter (will be removed)"
echo

#remove fastq files with 0 reads added to them:
for f in ${cutadapt_outputs}/plate??_*_${reads_name}_cutadapt_porechop.fastq; do
    [ -e "$f" ] || continue

    # NEW: strip out flagged multi-adapter reads before the empty check below
    if [ -s "$multi_ids" ]; then
        filtered="${f%.fastq}_filtered.fastq"
        seqkit grep -v -f "$multi_ids" "$f" -o "$filtered"
        mv "$filtered" "$f"
    fi

    lines=$(wc -l < "$f")
    if [ "$lines" -lt 4 ]; then
        rm -f "$f"
    fi
done

#Use find
find "${cutadapt_outputs}" -type f -name 'plate??_*.fastq' | while read -r plate_file_path; do
        echo "entered well cutadapt loop, plate_file_path = ${plate_file_path} "

        echo

        #Move the plate_ demultiplexed files we just made into directories based off the file names:
        plate_file_name="$(basename "$plate_file_path")"
        plate="${plate_file_name%%_*}"
        plate_dir="${cutadapt_outputs}/${plate}";
        mkdir -p "${plate_dir}"
        #mkdir -p "${cutadapt_outputs}/${plate}";
        mv "${plate_file_path}" "${plate_dir}/";
        #mv "${plate_file_path}" "${cutadapt_outputs}/${plate}/";
        echo "File sorted into directory 'plate_dir' : ${plate_dir}  "
        echo

        #Ok we want to execute the well-demultiplexing step once for each plate file. so include it in this loop:

        echo "we just demuxed cutadapt_outputs by PLATE, and made a new cutadapt_outputs/plate?? directory for each plate barcode. We are about to start demultiplexing by well with cutadapt.\n"

        echo "lets check that ${plate_dir} which we just set, exist and contains the correct number of items with ls"
        ls -1 "$plate_dir" | wc -l

        #DEMULTIPLEXING - STEP 3 - WELL
        #07-30-2026 CM update: change to O=18 (75%)
        #07-30-2026 CM update: change the 5' adapter search to cut from the RIGHTMOST match if multiple matches to the SAME 5' adapter are found in a single read
        echo "executing demultiplaexing step 3: cutadapt search for well barcodes! input file=${plate_dir}/${plate_file_name} "
        echo

        #70-31-2026 CM update: use an array to expand all the barcode names in well barcodes file because file: and ;rightmost arent compatible in cutadapt
        # NEW: work around a cutadapt bug (reproduced on 5.2 and current dev build) where
        # file:...;rightmost together crash with "unexpected keyword argument 'rightmost'".
        # Build individual -g "name=SEQ;rightmost" args from the fasta instead of using file:.
        declare -a well_g_args=()
        while read -r wb_name && read -r wb_seq; do
                wb_name="${wb_name#>}"
                well_g_args+=(-g "${wb_name}=${wb_seq};rightmost")
        done < "${well_barcodes}"

        cutadapt "${well_g_args[@]}" -O 18 --revcomp -e 0.15 --cores=0 --info-file ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_dir}/${plate_file_name} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log

        # END declare/rightmost workaround code


        # USE THIS ONE if you abandon the rightmost workaround:
        #cutadapt -g "file:${well_barcodes};rightmost" -O 18 -e 0.15 --cores=0 --info-file ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_dir}/${plate_file_name} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log;

        #cutadapt -g "GCGAGTCTTGT";rightmost -O 6 -e 0.15 --cores=0 --info-file ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_dir}/${plate_file_name} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log;
        #cutadapt -g file:${well_barcodes} -O 14 --action=lowercase --revcomp -e 0.15 -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_file_path} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log;

        # Remove empty well FASTQs
        for f in "${plate_dir}"/plate??_well??_*_${reads_name}_cutadapt_porechop.fastq; do
                [ -e "$f" ] || continue
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
                #[ -e "$plate_well_file_path" ] || continue
                # Extract the filename
                plate_well_file_name="$(basename "$plate_well_file_path")";

                # Extract the well ID (e.g., well01)
                well="${plate_well_file_name#*_}";
                well="${well%%_*}";
                #echo "current well = ${well}";

                # Create the well subdirectory
                plate_well_dir="$plate_dir/$well";
                mkdir -p "$plate_well_dir";

                # Move the file into the well subdirectory
                #echo "moving ${plate_well_file_path} to ${plate_well_dir}/"
                mv "${plate_well_file_path}" "${plate_well_dir}/";

                #while still looping thru values of plate and well, do cutadapt search for primers

                #DEMULTIPLEXING - STEP 4 - SEGMENT
                #okay demultiplex by plaque, input = plate-demuxed files; -O is smaller bc the primers are shorter
                cutadapt -a small=CTTTCGTACAACCGAGTAGG...CTCCTGAAGTATCTCACGCC -a medium=CGCTACGGCGGTATTGTC...GCTCACCAAGTAAGGTGTAGTAT -a large=TCGATGTTCAACTACTACGC...GCGAGACTCGCTTTGC -O 10 --revcomp -e 0.15 --cores=0 --info-file ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_well_dir}/${plate}_${well}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_well_dir}/${plate_well_file_name} > ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop.log;
                #cutadapt -a small=CTTTCGTACAACCGAGTAGG...CTCCTGAAGTATCTCACGCC -a medium=CGCTACGGCGGTATTGTC...GCTCACCAAGTAAGGTGTAGTAT -a large=TCGATGTTCAACTACTACGC...GCGAGACTCGCTTTGC -O 10 --action=lowercase --revcomp -e 0.15 -o ${plate_well_dir}/${plate}_${well}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_well_dir}/${plate_well_file_name} > ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop.log;

                #cutadapt -g file:${well_barcodes} -O 10 --revcomp -e 0.15 -o ${plate_well_dir}/${plate}_${well}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_well_dir}/${plate_well_file_name} > ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop.log;

                 #get read count in each file (before they are removed!):
                for fastq in ${plate_well_dir}/${plate}_${well}_*_${reads_name}_cutadapt_porechop.fastq;
                        do count=$( wc -l ${fastq} | awk '{print $1 / 4}');
                        echo "${cross},${count}" >> ${working_dir}/file_counts.csv;
                done

                # Remove empty segment FASTQs
                for f in ${plate_well_dir}/${plate}_${well}_*_${reads_name}_cutadapt_porechop.fastq; do
                        [ -e "$f" ] || continue
                        lines=$(wc -l < "$f")
                        if [ "$lines" -lt 4 ]; then
                                rm -f "$f"
                        fi
                done

                # If directory is now empty, delete it
                rmdir "${plate_well_dir}" 2>/dev/null

        done
done

echo "Dont forget!! I moved your previous cutadapt_outputs and porechop_outputs directories to cutadapt_outputs_from_${ts_cut} & porechop_outputs_from_${ts_pore} :)"
