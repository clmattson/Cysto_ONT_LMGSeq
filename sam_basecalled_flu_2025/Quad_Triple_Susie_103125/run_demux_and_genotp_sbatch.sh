#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=20:00:00
#SBATCH --mem=60G
#SBATCH --partition=low
#SBATCH --output=run_pipeline_%j.log
#SBATCH --error=run_pipeline_%j.err

echo "Running on $(hostname)"
echo "Job started at: $(date)"

###############################################
# Hybrid Input System:
# Priority: flags > prompts
# (Slurm does NOT export flag values, so we do
#  NOT use SBATCH env overrides.)
###############################################

# 1. Parse command-line flags
while getopts "d:r:p:w:c:l:s:e:" opt; do
    case $opt in
        d) DEMUXED_PATH="$OPTARG" ;;
        r) READS_PATH="$OPTARG" ;;
        p) PLATE_BARCODES="$OPTARG" ;;
        w) WELL_BARCODES="$OPTARG" ;;
        c) DO_PORECHOP="$OPTARG" ;;   # preserves "0" correctly
        l) MIN_LENGTH="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        e) LIB_PREP_DATABASE="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

echo "----------------------------------------"
echo " Welcome to the Demux + Genotyping Pipeline"
echo "----------------------------------------"

###############################################
# 2. Interactive prompts for missing values
# NOTE: we use:  [[ -z "${VAR+x}" ]]
# This checks if VAR is UNSET, not if it is "0".
###############################################

if [[ -z "${DEMUXED_PATH+x}" ]]; then
    while true; do
        read -p "Enter working/demuxed directory: " DEMUXED_PATH
        [[ -d "$DEMUXED_PATH" ]] && break
        echo "Directory not found. Try again."
    done
fi

if [[ -z "${READS_PATH+x}" ]]; then
    while true; do
        read -p "Enter path to raw reads: " READS_PATH
        [[ -f "$READS_PATH" ]] && break
        echo "Reads path not found. Try again."
    done
fi

if [[ -z "${PLATE_BARCODES+x}" ]]; then
    while true; do
        read -p "Enter plate barcode FASTA: " PLATE_BARCODES
        [[ -f "$PLATE_BARCODES" ]] && break
        echo "File not found. Try again."
    done
fi

if [[ -z "${WELL_BARCODES+x}" ]]; then
    while true; do
        read -p "Enter well barcode FASTA: " WELL_BARCODES
        [[ -f "$WELL_BARCODES" ]] && break
        echo "File not found. Try again."
    done
fi

if [[ -z "${DO_PORECHOP+x}" ]]; then
    while true; do
        read -p "Run porechop? (y/n): " ans
        case $ans in
            y|Y) DO_PORECHOP=1; break ;;
            n|N) DO_PORECHOP=0; break ;;
            *) echo "Please enter y or n." ;;
        esac
    done
fi

# NOTE: the previous version of this wrapper never actually asked for or passed a
# minimum-length value to demultiplexing.sh at all - it was silently missing from the
# call below, which means --minimum-length was being passed an empty value. Added here.
if [[ -z "${MIN_LENGTH+x}" ]]; then
    read -p "Enter minimum read length to keep [100]: " MIN_LENGTH
    MIN_LENGTH="${MIN_LENGTH:-100}"
fi

if [[ -z "${SAMPLE_LIST+x}" ]]; then
    while true; do
        read -p "Enter sample list CSV: " SAMPLE_LIST
        [[ -f "$SAMPLE_LIST" ]] && break
        echo "CSV file not found. Try again."
    done
fi

if [[ -z "${LIB_PREP_DATABASE+x}" ]]; then
    while true; do
        read -p "Enter library prep reference database FASTA: " LIB_PREP_DATABASE
        [[ -f "$LIB_PREP_DATABASE" ]] && break
        echo "File not found. Try again."
    done
fi

###############################################
# 3. Create timestamped pipeline log
###############################################
LOGFILE="${DEMUXED_PATH}/pipeline_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "----------------------------------------"
echo "Inputs confirmed:"
echo "  Working dir:       $DEMUXED_PATH"
echo "  Reads:             $READS_PATH"
echo "  Plate barcodes:    $PLATE_BARCODES"
echo "  Well barcodes:     $WELL_BARCODES"
echo "  Porechop:          $DO_PORECHOP"
echo "  Min length:        $MIN_LENGTH"
echo "  Sample list CSV:   $SAMPLE_LIST"
echo "  Lib prep database: $LIB_PREP_DATABASE"
echo "  Log file:          $LOGFILE"
echo "----------------------------------------"
echo "Note: demultiplexing.sh and genotyping.sh also each write their own stage-level"
echo "logs inside cutadapt_outputs/ and genotyping_outputs/ respectively - this file"
echo "is the overall pipeline-level record, not a replacement for those."
echo "----------------------------------------"

###############################################
# Step 1: demultiplexing.sh
###############################################
echo "Activating conda environment: cutadapt"
source ~/.bashrc
conda activate cutadapt

echo "Starting Step 1: demultiplexing_log_and_check.sh at $(date)"
/usr/bin/time -v bash demultiplexing_log_and_check.sh \
    -d "$DEMUXED_PATH" \
    -r "$READS_PATH" \
    -p "$PLATE_BARCODES" \
    -w "$WELL_BARCODES" \
    -c "$DO_PORECHOP" \
    -l "$MIN_LENGTH"
echo "Step 1 complete at $(date)"

###############################################
# Step 2: genotyping.sh
###############################################
echo "Loading module: usearch"
module load usearch

echo "Starting Step 2: genotyping_log_and_check.sh at $(date)"
# NOTE: previous versions of this wrapper called genotyping with "-e" for the sample
# list, but genotyping.sh's actual flag for that is "-s" (see genotyping.sh's own
# getopts string: d:s:l:). Fixed below - was silently mismatched before.
/usr/bin/time -v bash genotyping_log_and_check.sh \
    -d "$DEMUXED_PATH" \
    -s "$SAMPLE_LIST" \
    -l "$LIB_PREP_DATABASE"
echo "Step 2 complete at $(date)"

###############################################
# Step 3: Cutadapt log parser - PARKED FOR NOW
###############################################
# This step (double_parse_cutadapt_log_fast.python / batch_double_parse.sh) hasn't
# been updated to match the current demultiplexing.sh/genotyping.sh output structure
# yet. Left here, commented out, so it's easy to re-enable once that script is
# brought up to date - just uncomment and confirm the paths still match.
#
# CUTADAPT_DIR="${DEMUXED_PATH}/cutadapt_outputs"
#
# if [[ ! -d "$CUTADAPT_DIR" ]]; then
#     echo "ERROR: cutadapt_outputs not found at: $CUTADAPT_DIR"
#     exit 1
# fi
#
# cp double_parse_cutadapt_log_fast.python "$CUTADAPT_DIR/"
# cp batch_double_parse.sh "$CUTADAPT_DIR/"
# cd "$CUTADAPT_DIR" || exit 1
#
# echo "Starting Step 3: Cutadapt log parser wrapper at $(date)"
# /usr/bin/time -v bash batch_double_parse.sh
# echo "Step 3 complete at $(date)"

echo "Pipeline finished successfully at $(date)"

