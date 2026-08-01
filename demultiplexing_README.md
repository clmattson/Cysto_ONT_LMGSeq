README: Dual-Barcoded Nanopore Read Demultiplexing Pipeline
Overview

This script performs hierarchical demultiplexing and filtering of Nanopore reads generated from dual-PCR-barcoded libraries.

The pipeline:

Optionally runs Porechop to split chimeric reads.
Uses Cutadapt to identify and demultiplex:
Plate barcodes
Well barcodes
Segment-specific primer pairs
Filters reads below a user-specified minimum length.
Identifies and removes reads containing multiple distinct plate barcodes.
Organizes output FASTQ files into a nested plate/well directory structure.

The final output is a directory tree containing reads grouped by:

Plain Text
1
Plate
2
└── Well
3
└── Segment
Show more lines

This script is intended for Oxford Nanopore FASTQ data generated from the dual-barcoding workflow used within our group.

Before You Start
Important Assumptions

The script assumes:

Basecalling has already been completed.
Input reads are contained in a single FASTQ file.
Dorado trimming was disabled during basecalling.
Plate and well barcode FASTA files have already been prepared.
You are running on the HPC cluster rather than a login node.
Running on the HPC
Request Interactive Resources

This pipeline can be computationally intensive, particularly during:

Porechop processing
Cutadapt demultiplexing
Large FASTQ handling

Do not run this directly on a login node.

Start an interactive job:

Shell
1
srun -c 15 -t 40:00:00 --mem=50000 --partition=low --pty /bin/bash
Show more lines
Resource Explanation
Option	Meaning-c 15	Request 15 CPU cores
-t 40:00:00	Request 40 hours runtime
--mem=50000	Request 50 GB RAM
--partition=low	Submit to the low-priority partition
--pty /bin/bash	Open an interactive shell

For most sequencing runs these resources should be sufficient.

If processing unusually large datasets, increasing memory may be necessary.

Alternative: Batch Submission

For long unattended runs, consider submitting through sbatch.

Example:

Shell
1
#!/bin/bash
2
#SBATCH -c 15
3
#SBATCH --mem=50G
4
#SBATCH -t 40:00:00
5
#SBATCH -p low
6
 
7
bash demultiplexing_read_rename_filter_multiple_barcodes.sh \
8
-d /path/to/working_directory \
9
-r /path/to/reads.fastq \
10
-p /path/to/plate_barcodes.fasta \
11
-w /path/to/well_barcodes.fasta \
12
-c 1 \
13
-l 100
Show more lines

Submit with:

Shell
1
sbatch run_demux.sh
Show more lines
Conda Environment

This script requires software not available in a standard shell environment.

Activate the appropriate Conda environment before running:

Shell
1
conda activate cutadapt
Show more lines

Depending on future updates to the environment, the following programs should be available:

Shell
1
cutadapt
2
porechop
3
seqkit
4
``
Show more lines

Verify installation:

Shell
1
which cutadapt
2
which porechop
3
which seqkit
Show more lines

or

Shell
1
cutadapt --version
2
porechop --help
3
seqkit version
Show more lines

If any command is not found, the Conda environment was not activated correctly.

Input Files
Reads FASTQ

A single FASTQ containing all basecalled reads.

Example:

Plain Text
1
/group/project/run1/sup_basecall/FBC73506_fastq_pass.fastq
Show more lines
Plate Barcode FASTA

FASTA file containing plate barcode sequences.

Example:

Plain Text
fasta isn’t fully supported. Syntax highlighting is based on Plain Text.
1
>plate01
2
ACTGACTGACTG
3
 
4
>plate02
5
TGCATGCATGCA
Show more lines
Well Barcode FASTA

FASTA file containing well barcode sequences.

Example:

Plain Text
fasta isn’t fully supported. Syntax highlighting is based on Plain Text.
1
>well01
2
AACTTGGT
3
 
4
>well02
5
TTGGCCAA
Show more lines
Required Paths

All supplied paths should be absolute paths.

What is an absolute path?

An absolute path starts from the filesystem root (/) and uniquely identifies a file or directory.

Example:

Shell
1
/group/sldmunozgrp/project/reads.fastq
Show more lines

Not recommended:

Shell
1
../reads.fastq
Show more lines

Using absolute paths reduces mistakes and makes runs easier to reproduce.

Script Usage
General Command
Shell
1
bash demultiplexing_read_rename_filter_multiple_barcodes.sh \
2
-d <working_directory> \
3
-r <reads.fastq> \
4
-p <plate_barcodes.fasta> \
5
-w <well_barcodes.fasta> \
6
-c <0_or_1> \
7
-l <minimum_length>
Show more lines
Example Command
Shell
1
bash demultiplexing_read_rename_filter_multiple_barcodes.sh \
2
-d /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall \
3
-r /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall/FBC73506_fastq_pass_9ee158db_6bd36f99_0.fastq \
4
-p /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall/plate_barcodes.fasta \
5
-w /group/sldmunozgrp/cysto_LMGSeq08-25/Feb_cysto_flu/sup_basecall/plate_barcodes.fasta \
6
-c 0 \
7
-l 100
Show more lines
Command Line Flags
-d

Working directory.

This directory will contain:

Plain Text
1
cutadapt_outputs/
2
porechop_outputs/
3
file_counts.csv
4
 
Show more lines

Example:

Shell
1
-d /group/project/run1
Show more lines
-r

Input FASTQ file.

Example:

Shell
1
-r /group/project/run1/basecalled_reads.fastq
Show more lines
-p

Plate barcode FASTA.

Example:

Shell
1
-p /group/project/barcodes/plate_barcodes.fasta
Show more lines
-w

Well barcode FASTA.

Example:

Shell
1
-w /group/project/barcodes/well_barcodes.fasta
Show more lines
-c

Controls whether Porechop is executed.

-c 1

Run a fresh Porechop analysis.

Use when:

Starting from new reads.
No existing porechop outputs exist.
You want to regenerate chimeric-read splitting.

Example:

Shell
1
-c 1
Show more lines
-c 0

Skip Porechop and reuse an existing:

Plain Text
1
porechop_outputs/
Show more lines

directory.

Use when:

Testing parameter changes.
Re-running downstream demultiplexing.
Debugging.

Example:

Shell
1
-c 0
Show more lines

Important: If you specify -c 0 and porechop_outputs/ does not exist, the script will terminate with an error.

-l

Minimum read length after trimming.

Reads shorter than this value are discarded during Cutadapt processing.

Example:

Shell
1
-l 100
Show more lines
What Happens During the Run?
Step 1: Existing Outputs Are Archived

If present:

Plain Text
1
cutadapt_outputs/
2
porechop_outputs/
Show more lines

are renamed using timestamps.

Example:

Plain Text
1
cutadapt_outputs_from_07-30-2026_14.22
2
``
Show more lines

This prevents overwriting previous analyses.

Step 2: Porechop Splits Chimeric Reads (Optional)

Only performed when:

Shell
1
-c 1
Show more lines

The script:

Splits the FASTQ into chunks.
Runs Porechop on each chunk.
Collects statistics.
Fixes duplicate read IDs generated during splitting.
Recombines all processed reads into one FASTQ.

Output:

Plain Text
1
porechop_outputs/
Show more lines

Main file:

Plain Text
1
<reads_name>_porechop.fastq
Show more lines
Step 3: Plate Demultiplexing

Cutadapt identifies plate barcodes.

Output files:

Plain Text
1
plate01_*.fastq
2
plate02_*.fastq
3
...
Show more lines

Additional files:

Plain Text
1
plate_*_INFO.tsv
2
plate_*_cutadapt_porechop.log
Show more lines
Step 4: Multi-Adapter QC

The script examines the Cutadapt INFO file.

Reads containing more than one distinct plate barcode are flagged and removed.

Output:

Plain Text
1
<reads_name>_multi_adapter_ids.txt
Show more lines

This file contains the read IDs removed during QC.

Step 5: Well Demultiplexing

Within each plate:

Plain Text
1
plate01/
2
plate02/
3
...
Show more lines

the script searches for well barcodes.

New FASTQs are created and automatically sorted into well-specific subdirectories.

Example:

Plain Text
1
plate01/
2
├── well01/
3
├── well02/
4
├── well03/
5
└── ...
Show more lines
Step 6: Segment Identification

Within each well, Cutadapt searches for segment-specific primer pairs.

Categories:

Plain Text
1
small
2
medium
3
large
Show more lines

Output FASTQs are generated for each segment.

Example:

Plain Text
1
plate01/
2
└── well03/
3
├── plate01_well03_small.fastq
4
├── plate01_well03_medium.fastq
5
└── plate01_well03_large.fastq
Show more lines
Output Structure

A typical output tree may look like:

Plain Text
1
cutadapt_outputs/
2
├── plate01/
3
│ ├── plate01_well_INFO.tsv
4
│ ├── plate01_well.log
5
│ │
6
│ ├── well01/
7
│ │ ├── plate01_well01_small.fastq
8
│ │ ├── plate01_well01_medium.fastq
9
│ │ ├── plate01_well01_large.fastq
10
│ │ └── segment log and INFO files
11
│ │
12
│ └── well02/
13
│ └── ...
14
│
15
├── plate02/
16
│ └── ...
17
│
18
└── ...
Show more lines
Important Output Files
file_counts.csv

Contains read counts for generated segment FASTQs.

Useful for quick abundance summaries and QC.

Log Files

Generated at multiple stages:

Plain Text
1
*_cutadapt_porechop.log
Show more lines

and

Plain Text
1
*_porechop_chunk*.log
Show more lines

These contain:

adapter detection statistics
demultiplexing summaries
trimming statistics
troubleshooting information

Check these first if results appear unexpected.

INFO Files

Generated by Cutadapt:

Plain Text
1
*_INFO.tsv
Show more lines

These contain per-read adapter match information and are useful for:

troubleshooting barcode assignments
identifying reads with multiple adapter matches
validating barcode performance
Common Problems
"command not found"

Usually means:

Shell
1
conda activate cutadapt
Show more lines

was not run.

Porechop disabled but no outputs exist

If running:

Shell
1
-c 0
Show more lines

the directory:

Plain Text
1
porechop_outputs/
Show more lines

must already exist.

Otherwise the script will exit.

No Reads Assigned

Potential causes:

Incorrect barcode FASTA file.
Wrong barcode orientation.
Length filter too stringent.
Unexpected barcode sequence variation.

Review:

Plain Text
1
*_INFO.tsv
2
*_cutadapt_porechop.log
Show more lines

for adapter matching statistics.

Trailing Slashes in Paths

Avoid adding unnecessary trailing slashes:

Shell
1
/group/project/run1
Show more lines

preferred over:

Shell
1
/ group/project/run1/
Show more lines

The script generally tolerates them, but using consistent paths makes troubleshooting easier.

Recommended Workflow
Shell
1
srun -c 15 -t 40:00:00 --mem=50000 --partition=low --pty /bin/bash
2
 
3
conda activate cutadapt
4
 
5
bash demultiplexing_read_rename_filter_multiple_barcodes.sh \
6
-d /path/to/run_directory \
7
-r /path/to/reads.fastq \
8
-p /path/to/plate_barcodes.fasta \
9
-w /path/to/well_barcodes.fasta \
10
-c 1 \
11
-l 100
Show more lines

After completion, navigate to:

Shell
1
cd /path/to/run_directory/cutadapt_outputs
2
 
Show more lines

and explore the plate → well → segment hierarchy for your demultiplexed reads.
