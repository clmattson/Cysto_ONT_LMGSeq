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

#d path to the desired working dir
#r path to reads
#c plate_barcode.fasta file with path
#c well_barcode.fasta file with path
#s S genotyping locus reference path
#m M genotyping locus reference path
#l L genotyping locus reference path

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
cutadapt -a file:${plate_barcodes} -O 14 --revcomp -e 0.15 --minimum-length ${min_length} --info-file ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop_INFO.tsv -o ${cutadapt_outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${porechop_outputs}/${reads_name}_porechop.fastq > ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop.log
#cutadapt -a file:${plate_barcodes} -O 14 --action=lowercase --revcomp -e 0.15 -o ${cutadapt_outputs}/{name}_${reads_name}_cutadapt_porechop.fastq ${working_dir}/${reads_name}_porechop.fastq > ${cutadapt_outputs}/plate_${reads_name}_cutadapt_porechop.log

#remove fastq files with 0 reads added to them:
for f in ${cutadapt_outputs}/plate??_*_${reads_name}_cutadapt_porechop.fastq; do
    [ -e "$f" ] || continue
    lines=$(wc -l < "$f")
        if [ "$lines" -lt 4 ]; then
                rm -f "$f"
    fi
done

echo "completed demultiplexing step 2 - cutadapt plate identification!"
echo

#Use find
find "${cutadapt_outputs}" -type f -name 'plate??_*.fastq' | while read -r plate_file_path; do
        echo "entered cutadapt loop, plate_file_path = ${plate_file_path} "

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
        cutadapt -g file:${well_barcodes};rightmost -O 14 --revcomp -e 0.15 --info-file ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_dir}/${plate}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_dir}/${plate_file_name} > ${plate_dir}/${plate}_well_${reads_name}_cutadapt_porechop.log;
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
        for plate_well_file_path in "${plate_dir}"/plate??_well??_*.fastq; do
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
                cutadapt -a small=CTTTCGTACAACCGAGTAGG...CTCCTGAAGTATCTCACGCC -a medium=CGCTACGGCGGTATTGTC...GCTCACCAAGTAAGGTGTAGTAT -a large=TCGATGTTCAACTACTACGC...GCGAGACTCGCTTTGC -O 10 --revcomp -e 0.15 --info-file ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop_INFO.tsv -o ${plate_well_dir}/${plate}_${well}_{name}_${reads_name}_cutadapt_porechop.fastq ${plate_well_dir}/${plate_well_file_name} > ${plate_well_dir}/${plate}_${well}_segment_${reads_name}_cutadapt_porechop.log;
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
