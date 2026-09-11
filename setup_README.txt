# Setup Instructions

This document walks through everything needed to get my demultiplexing and genotyping
pipeline running from scratch: raw sequencing data on Rosalind, all the way through to
genotyped output on the HPC.


## Step 1: Basecalling raw data with Dorado (on Rosalind)

Raw sequencing data comes off the MinION as either pod5 or fast5 files. The best thing to do is set MinKNOW to generate .pod5 files, but it is possible to start with fast5s instead. If you're starting with pod5s, you can skip straight to the basecalling command below.

### a) Convert fast5 to pod5 (skip this if you already have pod5s)

fast5 data is usually split into `fast5_pass` and `fast5_fail` folders by fast basecalling algorithm run during sequencing. Convert each separately since we will want to re-run a more accurate basecalling algorithm on the complete data:
```
pod5 convert fast5 --output converted_fail_pod5s.pod5 fast5_fail/*/*.fast5
pod5 convert fast5 --output converted_pass_pod5s.pod5 fast5_pass/*/*.fast5
mkdir converted_pod5s
mv converted_*pod5s.pod5 converted_pod5s
cd converted_pod5s
```
If you want, you can check the read count in each converted file matches what you'd expect with `pod5 inspect
reads <file>.pod5 | wc -l`.

### b) Merge into a single pod5 file

```
pod5 merge converted_pass_pod5s.pod5 converted_fail_pod5s --output <run_name>_all.pod5
```
You can again check the merged read count is the sum of the two inputs, if you like. You can also delete the intermediate pod5 files since they're large and no longer needed after merging.

### c) Download the `sup` basecalling model (only needed once per model)

Run this command in the new directory with your pod5 data:
```
~/dorado/bin/dorado download --model sup --data <run_name>_all.pod5
```

### d) Run basecalling

```
~/dorado/bin/dorado basecaller -v -r --emit-fastq --min-qscore 8 --disable-read-splitting \
  --no-trim --output-dir <output_dir> sup <pod5_folder>
```
What each flag does:
- `-v` - verbose logging
- `-r` - recursive, so it finds pod5 files in subfolders too
- `--emit-fastq` - output fastq instead of the default BAM. This is required - my pipeline
  needs a fastq to start, not a BAM.
- `--min-qscore 8` - only keep reads with an average quality score of 8 or higher (this is
  the standard ONT pass/fail quality cutoff)
- `--disable-read-splitting` - don't let dorado auto-split reads on suspected internal
  adapters. Read splitting is handled downstream by porechop instead.
- `--no-trim` - don't trim adapters during basecalling. Adapter/barcode trimming is handled
  downstream by porechop and cutadapt, so trimming here would interfere with that.
- `sup` - the model name (must match what you downloaded in step c)
- last argument - the folder containing your pod5 files

### e) Combine into one fastq

Dorado may write out more than one fastq file depending on version/settings. My pipeline
needs a single combined fastq of all quality-passed reads to start from:
```
cat <output_dir>/*.fastq > <run_name>_calls.fastq
```
That combined fastq is what you'll transfer to the HPC and use as the `-r` input to
`demultiplexing.sh`.

## Step 2: Software setup

### a) Create the conda environment

I've exported my environment to `ont-lmgseq.yml` using `conda env export --from-history`,
which only lists the packages I explicitly installed.

To recreate it run the following in your terminal:
```
conda env create -f ont-lmgseq.yml
conda activate ont-lmgseq
```
The environment has the following tools: cutadapt, seqkit, seqtk, blast, biopython, and porechop** (extra setup required for porechop - see below**).
usearch is separate - it's installed as an HPC module, not through conda. Load it with:
```
module load usearch
```

Both `demultiplexing.sh` and `genotyping.sh` check for their required tools on PATH before
doing anything else, so if you forget one of these steps you'll get an error message.

### b) Install porechop

Porechop needs to be installed from source, because it needs a custom adapter file (see
part c) swapped in before installing - a plain `pip install porechop` or `conda install
porechop` will give you the stock version with the default adapter sequences.

The following installation instructions are based on the porechop github: https://github.com/rrwick/porechop.

a) Important: First, make sure the `ont-lmgseq` conda environment is active before you install it:
```
conda activate ont-lmgseq
```
Porechop's executable gets installed directly into whichever environment is active at the time of installation. For this reason, the `ont-lmgseq` conda environment also has to be active every time you run porechop, even though it isn't a conda package. This installation method worked reliably for me but I will mention another option below that should work***, but I have not tested.

b) Fetch the porechop tool and all files by cloning the porechop repo to your desired location:
```
git clone https://github.com/rrwick/Porechop.git
cd Porechop
```

### c) Add the custom adapter sequences

The `adapters.py` file in the Porechop source defines what sequences porechop will look for and trim. Our version only needs the PCR2 landing pad sequences (forward and reverse-complement), not the default, which is a list of Nanopore kit adapters.

Copy the custom `adapters.py` from this repo over the one in the cloned source:
```
cp /path/to/repo/porechop_setup/adapters.py Porechop/porechop/adapters.py
```

Then continue with the porechop installation. Inside your cloned Porechop/ directory run:
```
python3 setup.py install
```

Confirm it worked:
```
porechop -h
python3 -c "from porechop.adapters import ADAPTERS; print([a.name for a in ADAPTERS])"
```
The second command should print `PCR2_landing_pads` and `PCR2_landing_pads_revcomp`, not a long list of Nanopore kit names. If you ever need to change the adapter sequences again, manually edit `adapters.py` and rerun `python3 setup.py install` from inside the Porechop folder with `ont-lmgseq` still active.

*** Alternative option: theoretically you should be able to use `pip3 install local/path/to/Porechop` AFTER you clone the repo and replace adapters.py, and would then avoid the need to use a conda environment, but I haven't tried it since the environment is needed for the other tools anyway.



## Step 3: Understanding the scripts

Each script prints its flags near the top as well (if you run the script on its own) and also have the info stored as comments near the top of the script files.

a) `demultiplexing.sh` - runs porechop then cutadapt to demultiplex reads by plate and well
   inputs:
    -d = absolute path to the desired working dir
    -r = absolute path to reads
    -p = <plate_barcode.fasta> file with plate barcodes fasta with absolute path
    -w = <well_barcode.fasta> file with well barcodes fasta with absolute path
    -c = do_porechop: 0 or 1 for whether or not to re-do the porechop step
    -l = min_length: minimum desired length cutoff for filtering reads

  outputs: 
    porechop_outputs
    data: 
    logfiles:
    cutadapt_outputs
    data:  
    logfiles:

b) `genotyping.sh` - runs usearch (and optionally blastn) to assign a segment and strain identity to each demultiplexed well
    
  inputs:
    -d - path to the working directory that CONTAINS cutadapt_outputs and porechop_outputs - this is the SAME working directory you passed to the demultiplexing script
    -s - sample list  - CSV(!!) file (wih path) with all samples: ie barcode, cross, parent 1, parent 2"
    -l - library prep database = a fasta file with one reference genome each for ALL strains included in the library prep. Indiv sample strain assignment datbases get built from this. Include one reference per strain per genome segment and name them like this: 
      >Strain_Seg 
      ATCGATCGATGC....
    -b - do_blast: 0 or 1 only - enter 1 if you want to blast the reads that usearch didn't match to any ref. 

  outputs:   
    genotyping_outputs
      coinfection/positive/negative/misassigned folders with .b6 output files from usearch, sorted by experiment type
        data:
        log files:
      strain_assignment_output
        data:  
        log files:

c) `sbatch_demultiplexing.sh` / `sbatch_genotyping.sh` - sbatch wrappers for running the steps separately as a slurm job instead of interactively
  inputs and outputs are essentially the same as above, plus a log for the slurm job
d) `run_pipeline_sbatch.sh` - runs both steps back to back in a single sbatch 
  inputs and outputs are essentially the same as above but combined (so you must provide all of the flags for both), plus a log for the slurm job


## Step 4: Get your data and supporting files

You'll need all of the following before running the pipeline:

### a) Basecalled reads
A single fastq file containing all quality-passed reads (see Step 2e). It is best to name it something meaningful, like "Experiment1_MM_DD_YYYY_Name.fastq" if possible

### b) Barcode files
Two fasta files are alread privided: plate_barcodes.fasta and well_barcodes.fasta. These have just the sequences of just the unique identifier portion of the LMGSeq PCR primer sequences. Each fasta header is the barcode's name (e.g. `>plate01`, `>well01`) and the sequence below it is the actual barcode sequence.

### c) Sample list CSV
A CSV listing every plate/well combination and what it is: sample type (positive, negative, coinfection, or misassigned) and the parent strain(s) involved.

The sample list must have the following columns. `plate` and `well` refer to the barcode numbers assigned to the sample at library prep. `parentX` is the name of the parent REFERENCE in the database, and `parentX_label` can be another name if you would like. If a particular sample does not have a relevant entry for a particular column, please make sure leave it blank.

i) plate,well,type,parent1,parent2,parent3,parent4,parent1_label,parent2_label,parent3_label,parent4_label

ii) Rows 1-3 show the layout for coinfection samples with TWO parents:
plate01,well01,coinfection,HK68,PAN99,,,HK68,PAN99,,
plate01,well02,coinfection,HK68,PAN99,,,HK68,PAN99,,
plate01,well03,coinfection,HK68,PAN99,,,HK68,PAN99,,

iii) Rows 4-6 show the layout for coinfection samples with THREE parents:
plate03,well10,coinfection,PAN99,SI86,TX12,,PAN99,SI86,TX12,
plate03,well11,coinfection,PAN99,SI86,TX12,,PAN99,SI86,TX12,
plate03,well12,coinfection,PAN99,SI86,TX12,,PAN99,SI86,TX12,

iv) Rows 6-9 show the layout for coinfection samples with FOUR parents:
plate04,well70,coinfection,CH83PR8,SI86,TX12,CA09,CH83PR8,SI86,TX12,CA09
plate04,well71,coinfection,CH83PR8,SI86,TX12,CA09,CH83PR8,SI86,TX12,CA09
plate04,well72,coinfection,CH83PR8,SI86,TX12,CA09,CH83PR8,SI86,TX12,CA09

v) Rows 10 and 11 show the layout for positive control samples (ie not coinfections - single strain stocks, supernatants/lysates, or plaques):
plate15,well94,positive,PAN99,,,,PAN99,,,,
plate15,well95,positive,TX12,,,,TX12,,,,

vi) Rows 12 and 13 show the layout for negative control samples (ie water or buffer controls):
plate15,well96,negative,,,,,NFW_1,,,,
plate25,well78,negative,,,,,NFW_2,,,,

vii) Rows 14-16 show coinfection samples where the parent and parent label differ, because I included coinfections with different mutant strains but wanted to usearch them all against the reference for phi6:
plate20,well01,coinfection,PHI6,CA68,4267LP1,CA68
plate21,well02,coinfection,PHI6,CA68,4267LP3,CA68
plate22,well03,coinfection,PHI6,CA68,4267LP5,CA68


So all of these rows would get combined into one file that looks like this:



### d) Library prep reference database
A single fasta file containing every reference genome used anywhere in this library prep.
One reference per strain per genome segment, named `>Strain_Segment` - for example:
```
>Strain_Seg 
TCGATCGATGCATCGATCGATGCATCGATCGATGCATCGATCGATGCA
>PAN99_M
CGATCGATGCATCGATCGATGCATCGATCGATGCATCGATCGATCGAT
>CA09_NP
ATCGATCGATGCATCGATCGATGCATCGATCGATGCATCGATCGATGC



Once you have all four of these files, you're ready to run `demultiplexing.sh` followed by
`genotyping.sh`
