#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH --partition=low
#SBATCH --output=demultiplexing_%j.log
#SBATCH --error=demultiplexing_%j.err

echo "Running on $(hostname)"
echo "Job started at: $(date)"

###############################################
# Hybrid Input System:
# Priority: flags > prompts
# (Slurm does NOT export flag values, so we do
#  NOT use SBATCH env overrides.)
###############################################

while getopts "d:r:p:w:c:l:" opt; do
    case $opt in
        d) DEMUXED_PATH="$OPTARG" ;;
        r) READS_PATH="$OPTARG" ;;
        p) PLATE_BARCODES="$OPTARG" ;;
        w) WELL_BARCODES="$OPTARG" ;;
        c) DO_PORECHOP="$OPTARG" ;;   # preserves "0" correctly
        l) MIN_LENGTH="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

echo "----------------------------------------"
echo " Demultiplexing step"
echo "----------------------------------------"

###############################################
# Interactive prompts for missing values
# NOTE: we use:  [[ -z "${VAR+x}" ]]
# This checks if VAR is UNSET, not if it is "0".
###############################################

if [[ -z "${DEMUXED_PATH+x}" ]]; then
    while true; do
        read -p "Enter working directory: " DEMUXED_PATH
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

if [[ -z "${MIN_LENGTH+x}" ]]; then
    read -p "Enter minimum read length to keep [100]: " MIN_LENGTH
    MIN_LENGTH="${MIN_LENGTH:-100}"
fi

echo "----------------------------------------"
echo "Inputs confirmed:"
echo "  Working dir:      $DEMUXED_PATH"
echo "  Reads:            $READS_PATH"
echo "  Plate barcodes:   $PLATE_BARCODES"
echo "  Well barcodes:    $WELL_BARCODES"
echo "  Porechop:         $DO_PORECHOP"
echo "  Min length:       $MIN_LENGTH"
echo "----------------------------------------"

echo "Activating conda environment: cutadapt"
source ~/.bashrc
conda activate cutadapt

echo "Starting demultiplexing_log2_and_check.sh at $(date)"
/usr/bin/time -v bash demultiplexing_log2_and_check.sh \
    -d "$DEMUXED_PATH" \
    -r "$READS_PATH" \
    -p "$PLATE_BARCODES" \
    -w "$WELL_BARCODES" \
    -c "$DO_PORECHOP" \
    -l "$MIN_LENGTH"
echo "demultiplexing.sh complete at $(date)"

