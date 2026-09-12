#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH --partition=low
#SBATCH --output=genotyping_%j.log
#SBATCH --error=genotyping_%j.err

echo "Running on $(hostname)"
echo "Job started at: $(date)"

###############################################
# Hybrid Input System:
# Priority: flags > prompts
###############################################

while getopts "d:s:l:b:" opt; do
    case $opt in
        d) DEMUXED_PATH="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        l) LIB_PREP_DATABASE="$OPTARG" ;;
	b) DO_BLAST="$OPTARG" ;;
	*) echo "Invalid option"; exit 1 ;;
    esac
done

echo "----------------------------------------"
echo " Genotyping step"
echo "----------------------------------------"

if [[ -z "${DEMUXED_PATH+x}" ]]; then
    while true; do
        read -p "Enter working directory (the same one used for demultiplexing.sh, containing cutadapt_outputs): " DEMUXED_PATH
        [[ -d "$DEMUXED_PATH" ]] && break
        echo "Directory not found. Try again."
    done
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

if [[ -z "${DO_BLAST+x}" ]]; then
    while true; do
        read -p "Run blast? (y/n): " ans
        case $ans in
            y|Y) DO_BLAST=1; break ;;
            n|N) DO_BLAST=0; break ;;
            *) echo "Please enter y or n." ;;
        esac
    done
fi

echo "----------------------------------------"
echo "Inputs confirmed:"
echo "  Working dir:       $DEMUXED_PATH"
echo "  Sample list CSV:   $SAMPLE_LIST"
echo "  Lib prep database: $LIB_PREP_DATABASE"
echo "  Do blast: $DO_BLAST"
echo "----------------------------------------"

echo "Loading usearch"
module load usearch

echo "Starting genotyping_log2_and_check.sh at $(date)"
/usr/bin/time -v bash genotyping_log2_and_check.sh \
    -d "$DEMUXED_PATH" \
    -s "$SAMPLE_LIST" \
    -l "$LIB_PREP_DATABASE"
    -b "$DO_BLAST"
echo "genotyping_log_and_check.sh complete at $(date)"

