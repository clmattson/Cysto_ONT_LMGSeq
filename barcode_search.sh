#!/bin/bash

#to do:


#get input:

#

#-p ${parent_dir} directory containing data and (#sup_basecall)
#-d sample_suffix - starts with an underscore, example: "_norm3_meta_spades"
#-s ${host_s} S segment reference file & path
#-m ${host_m} M segment reference file & path
#-l ${host_l} L segment reference file & path


#parent_path=''
sample_suffix=''
#host_s=''
#host_m=''
#host_l=''


print_usage() {
  printf "Usage: ..."
}

while getopts d: flag
do
    case "${flag}" in
d) sample_suffix=${OPTARG};;
#s) host_s=${OPTARG};;
#m) host_m=${OPTARG};;
#l) host_l=${OPTARG};;

    esac
done


#cd ${parent_path} || exit

echo "working shell directory:"
pwd
echo




#okay % classified looking sad for the plaque primers. lets get some code that searches for the barcodes and identifies the % of reads containing an exact match

file="dorado_barcode_arrs/plaque_barcodes.fasta"


#barcode_matches_by_plaque=()
barcode_S_matches_by_plaque=()
barcode_M_matches_by_plaque=()
barcode_L_matches_by_plaque=()
plaque_barcodes_fasta="dorado_barcode_arrs/plaque_barcodes.fasta"
reads_fastq="fastq_pass/FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq"

total_reads=$( wc -l FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq | awk '{print $1/4}')

for plaque_barcode_seq in $(grep -v "^>" "$plaque_barcodes_fasta" ); do
  #get line in barcodes fasta with sequwnce
  seq_line=$(sed -n "/${plaque_barcode_seq}/=" ${plaque_barcodes_fasta});   
  echo "seq_line = ${seq_line}";
  #get line in barcodes fasta with header name
  header_line=$((seq_line -1));  
  #create variables for the barcoode sequence + the GSP fwd primer
  plaque_S_seq=$("${plaque_barcode_seq}CTTTCGTACAACCGAGTAGG");
  plaque_M_seq=$("${plaque_barcode_seq}CGCTACGGCGGTATTGTC");
  plaque_L_seq=$("${plaque_barcode_seq}TCGATGTTCAACTACTACGC");
  #search reads files for the above combos
  reads_with_S_barcode=$(grep -c "${plaque_S_seq}" ${reads_fastq});
  reads_with_M_barcode=$(grep -c "${plaque_M_seq}" ${reads_fastq});
  reads_with_L_barcode=$(grep -c "${plaque_L_seq}" ${reads_fastq});
  #print them out
  echo "reads_with_S_barcode = ${reads_with_S_barcode}";
  echo "reads_with_M_barcode = ${reads_with_M_barcode}";
  echo "reads_with_L_barcode = ${reads_with_L_barcode}";
  #total_reads=$( wc -l FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq | awk '{print $1/4}');
  echo "total_reads = ${total_reads}";
  #calculate % of total reads matvhed to the currect barcode, split by segment too
  percent_S_reads_with_barcode=$(echo "scale=5 ; $reads_with_S_barcode/$total_reads"| bc);
  percent_M_reads_with_barcode=$(echo "scale=5 ; $reads_with_M_barcode/$total_reads"| bc);
  percent_L_reads_with_barcode=$(echo "scale=5 ; $reads_with_L_barcode/$total_reads"| bc);
  #get header name
  header=$(sed -n "${header_line}p" ${plaque_barcodes_fasta}); 
  #print them out
  echo "${percent_S_reads_with_barcode} of the reads contain an exact S segment match to the barcode for ${header}";
  echo "${percent_M_reads_with_barcode} of the reads contain an exact M segment match to the barcode for ${header}";
  echo "${percent_L_reads_with_barcode} of the reads contain an exact L segment match to the barcode for ${header}";
  #append the fraction to an array
  barcode_S_matches_by_plaque+=("${percent_S_reads_with_barcode}");
  barcode_M_matches_by_plaque+=("${percent_M_reads_with_barcode}");
  barcode_L_matches_by_plaque+=("${percent_L_reads_with_barcode}");
  #print the array
  echo "${barcode_S_matches_by_plaque[@]}";
  echo "${barcode_M_matches_by_plaque[@]}";
  echo "${barcode_L_matches_by_plaque[@]}";

  #get num reads for current barcode in fastq output by dorado demux:
  plaque_num="${header#*plaque}"

  dorado_M_read_count=$( wc -l demux_plaque_M/*plaque_M_barcode${plaque_num}.fastq | awk '{print $1/4}');
  dorado_S_read_count=$( wc -l demux_plaque_S/*plaque_S_barcode${plaque_num}.fastq | awk '{print $1/4}');
  dorado_L_read_count=$( wc -l demux_plaque_L/*plaque_L_barcode${plaque_num}.fastq | awk '{print $1/4}');
 
  echo;
  echo;
done




barcode_matches_by_plaque=()
barcode_S_matches_by_plaque=()
barcode_M_matches_by_plaque=()
barcode_L_matches_by_plaque=()
dorado_barcode_S_matches_by_plaque=()
dorado_barcode_M_matches_by_plaque=()
dorado_barcode_L_matches_by_plaque=()
dorado_cross_plaque_barcode_S_matches_by_plaque=()
dorado_cross_plaque_barcode_M_matches_by_plaque=()
dorado_cross_plaque_barcode_L_matches_by_plaque=()
barcode_S_matches_by_cross_plaque=()
barcode_M_matches_by_cross_plaque=()
barcode_L_matches_by_cross_plaque=()
plaque_barcodes_fasta="dorado_barcode_arrs/plaque_barcodes.fasta"
reads_fastq="fastq_pass/FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq"
cross_barcodes_fasta="dorado_barcode_arrs/cross_barcodes.fasta"

total_reads=$( wc -l FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq | awk '{print $1/4}')





#loops thru lines of barcode file that dont contain (-V) the > symbiol, aka the sequences in the fasta file
for plaque_barcode_seq in $(grep -v "^>" "$plaque_barcodes_fasta" ); do
#get line in barcodes fasta with sequwnce
seq_line=$(sed -n "/${plaque_barcode_seq}/=" ${plaque_barcodes_fasta});      
echo "seq_line = ${seq_line}";  
#get line in barcodes fasta with fasta name
header_line=$((seq_line -1)); 
#search reads files for the above barcode seqs, spit out subset of reads containing the sequence:
grep -B 1 -A 2 --no-group-separator "${plaque_barcode_seq}" ${reads_fastq} > reads_with_barcode.fastq;
echo "grepped the full reads fastq on barcode sequence ${plaque_barcode_seq} and output matching reads to reads_with_barcode_tmp.fastq";
num_reads_with_barcode=$( wc -l reads_with_barcode.fastq | awk '{print $1/4}')
#num_reads_with_barcode=$(grep -c "${plaque_barcode_seq}" ${reads_fastq});   
echo "reads_with_barcode = ${num_reads_with_barcode}";
echo "total_reads = ${total_reads}";   
percent_reads_with_barcode=$(echo "scale=5 ; $num_reads_with_barcode/$total_reads"| bc);  
  #get header name
header=$(sed -n "${header_line}p" ${plaque_barcodes_fasta});    
echo "${percent_reads_with_barcode} fraaction of the reads contain an exact match to the barcode for ${header}";   
#append fraction to array
barcode_matches_by_plaque+=("${percent_reads_with_barcode}");   
echo "{$barcode_matches_by_plaque[@]}";   

#create variables for the barcoode sequence + the GSP fwd primer
plaque_S_primer="CTTTCGTACAACCGAGTAGG";
plaque_M_primer="CGCTACGGCGGTATTGTC";
plaque_L_primer="TCGATGTTCAACTACTACGC";
#search reads files for the above combos
num_reads_with_S_primer=$(grep -c "${plaque_S_primer}" reads_with_barcode.fastq);
num_reads_with_M_primer=$(grep -c "${plaque_M_primer}" reads_with_barcode.fastq);
num_reads_with_L_primer=$(grep -c "${plaque_L_primer}" reads_with_barcode.fastq);
#print them out
echo "reads_with_S_primer = ${num_reads_with_S_primer}";
echo "reads_with_M_primer = ${num_reads_with_M_primer}";
echo "reads_with_L_primer = ${num_reads_with_L_primer}";

#calculate % of total reads matvhed to the currect barcode, split by segment too
percent_S_reads_with_barcode=$(echo "scale=5 ; $num_reads_with_S_primer/$total_reads"| bc);
percent_M_reads_with_barcode=$(echo "scale=5 ; $num_reads_with_M_primer/$total_reads"| bc);
percent_L_reads_with_barcode=$(echo "scale=5 ; $num_reads_with_L_primer/$total_reads"| bc);
#print them out
echo "${percent_S_reads_with_barcode} of the reads contain an exact S segment match to the barcode for ${header}";
echo "${percent_M_reads_with_barcode} of the reads contain an exact M segment match to the barcode for ${header}";
echo "${percent_L_reads_with_barcode} of the reads contain an exact L segment match to the barcode for ${header}";
#append the fraction to an array
barcode_S_matches_by_plaque+=("${num_reads_with_S_primer}");
barcode_M_matches_by_plaque+=("${num_reads_with_M_primer}");
barcode_L_matches_by_plaque+=("${num_reads_with_L_primer}");
#print the array
echo "barcode matches by plaque found by grep:"
echo "${barcode_S_matches_by_plaque[@]}";
echo "${barcode_M_matches_by_plaque[@]}";
echo "${barcode_L_matches_by_plaque[@]}";

#get num reads for current barcode in fastq output by dorado demux:
plaque_num="${header#*plaque}"

dorado_S_read_count=$( wc -l demux_plaque_S/*plaque_S_barcode${plaque_num}.fastq | awk '{print $1/4}');
dorado_M_read_count=$( wc -l demux_plaque_M/*plaque_M_barcode${plaque_num}.fastq | awk '{print $1/4}');
dorado_L_read_count=$( wc -l demux_plaque_L/*plaque_L_barcode${plaque_num}.fastq | awk '{print $1/4}');

dorado_barcode_S_matches_by_plaque+=("${dorado_S_read_count}");
dorado_barcode_M_matches_by_plaque+=("${dorado_M_read_count}");
dorado_barcode_L_matches_by_plaque+=("${dorado_L_read_count}");

echo "barcode matches by plaque found by dorado:"
echo "${dorado_barcode_S_matches_by_plaque[@]}";
echo "${dorado_barcode_M_matches_by_plaque[@]}";
echo "${dorado_barcode_L_matches_by_plaque[@]}";

echo
echo "Dorado found ${dorado_S_read_count} reads with the ${header} barcode and the S primer | grep found ${num_reads_with_S_primer} reads with the S primer among ${num_reads_with_barcode} reads with the ${header} barcode."
echo "Dorado found ${dorado_M_read_count} reads with the ${header} barcode and the M primer | grep found ${num_reads_with_M_primer} reads with the M primer among ${num_reads_with_barcode} reads with the ${header} barcode."
echo "Dorado found ${dorado_L_read_count} reads with the ${header} barcode and the L primer | grep found ${num_reads_with_L_primer} reads with the L primer among ${num_reads_with_barcode} reads with the ${header} barcode."
echo;   
echo; 
done |& tee barcoding_loop_log.out



  for plaque_barcode_seq in $(grep -v "^>" "$plaque_barcodes_fasta" ); do   seq_line=$(sed -n "/${plaque_barcode_seq}/=" ${plaque_barcodes_fasta});      echo "seq_line = ${seq_line}";   header_line=$((seq_line -1));      reads_with_barcode=$(grep -c "${plaque_barcode_seq}" ${reads_fastq});   echo "reads_with_barcode = ${reads_with_barcode}";   echo "total_reads = ${total_reads}";   percent_reads_with_barcode=$(echo "scale=5 ; $reads_with_barcode/$total_reads"| bc);   header=$(sed -n "${header_line}p" ${plaque_barcodes_fasta});    echo "${percent_reads_with_barcode}% of the reads contain an exact match to the barcode for ${header}";   barcode_matches_by_plaque+=("${percent_reads_with_barcode}"); echo "${barcode_matches_by_plaque[@]}";   echo;   echo; done



#process the fully demuxed reads - cross and plaque


 declare -a arr=("S" "M" "L")

#run in sup_basecall
for i in "${arr[@]}"; do 
  for j in demux_plaque_${i}/demux_cross_plaque_${i}_barcode*; do 
    for reads_fastq in ${j}/*.fastq; do 
      echo; 
      echo "j = ${j}"; 
      reads_name="${reads_fastq##*/}"; 
      reads_name="${reads#*_}"; 
      echo "plaque_barcode_${j#*barcode}_${reads}"; 




####RERUN#############3

barcode_matches_by_plaque=()
barcode_S_matches_by_plaque=()
barcode_M_matches_by_plaque=()
barcode_L_matches_by_plaque=()
dorado_barcode_S_matches_by_plaque=()
dorado_barcode_M_matches_by_plaque=()
dorado_barcode_L_matches_by_plaque=()
dorado_cross_plaque_barcode_S_matches_by_plaque=()
dorado_cross_plaque_barcode_M_matches_by_plaque=()
dorado_cross_plaque_barcode_L_matches_by_plaque=()
barcode_S_matches_by_cross_plaque=()
barcode_M_matches_by_cross_plaque=()
barcode_L_matches_by_cross_plaque=()
plaque_barcodes_fasta="dorado_barcode_arrs/plaque_barcodes.fasta"
reads_fastq="fastq_pass/FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq"
cross_barcodes_fasta="dorado_barcode_arrs/cross_barcodes.fasta"
echo "plaque_num,cross_num,plaque_matches_all,s_cross_plaque_dorado,s_cross_plaque_grep,m_cross_plaque_dorado,m_cross_plaque_grep,l_cross_plaque_dorado,l_cross_plaque_grep" > dorado_vs_grep_barcode_matches.csv

total_reads=$( wc -l fastq_pass/FBC73506_fastq_pass_d5fa85e0_b727e37a_0.fastq | awk '{print $1/4}')

#echo "plaque_num,cross_num,s_plaque_dorado,s_plaque_grep,m_plaque_dorado,m_plaque_grep,l_plaque_dorado,l_plaque_grep,s_cross_plaque_dorado,s_cross_plaque_grep,m_cross_plaque_dorado,m_cross_plaque_grep,l_cross_plaque_dorado,l_cross_plaque_grep" > dorado_vs_grep_barcode_matches.csv

#################################################################
      

#get fasta lines without the header symbol >, loop thru all of them:
for plaque_barcode_seq in $(grep -v "^>" "$plaque_barcodes_fasta" ); do
  #get line in barcodes fasta with sequence
  plaque_seq_line=$(sed -n "/${plaque_barcode_seq}/=" ${plaque_barcodes_fasta});      
  echo "plaque_seq_line = ${plaque_seq_line}";  
  #get line in barcodes fasta with fasta name
  plaque_header_line=$((plaque_seq_line -1)); 
  #search reads FASTQ  for the above barcode seqs, -BA get lines above and below match; spit out subset of reads containing the sequence:
  grep -B 1 -A 2 --no-group-separator "${plaque_barcode_seq}" ${reads_fastq} > reads_with_plaque_barcode.fastq;
  echo "grepped the full reads fastq on barcode sequence ${plaque_barcode_seq} and output matching reads to reads_with_plaque_barcode_tmp.fastq";
  num_reads_with_plaque_barcode=$( wc -l reads_with_plaque_barcode.fastq | awk '{print $1/4}')
  #num_reads_with_plaque_barcode=$(grep -c "${plaque_barcode_seq}" ${reads_fastq});   
  echo "reads_with_plaque_barcode = ${num_reads_with_plaque_barcode}";
  echo "total_reads = ${total_reads}";   
  percent_reads_with_plaque_barcode=$(echo "scale=5 ; $num_reads_with_plaque_barcode/$total_reads"| bc);  
  #get header name
  plaque_header=$(sed -n "${plaque_header_line}p" ${plaque_barcodes_fasta});    
  echo "${percent_reads_with_plaque_barcode} fraaction of the reads contain an exact match to the barcode for ${plaque_header}";   
  #append fraction to array
  barcode_matches_by_plaque+=("${percent_reads_with_plaque_barcode}");   
  echo "{$barcode_matches_by_plaque[@]}";   

  #now lets search for cross barcode before we search for primers
  #repeat the loop and search steps
  #get fasta lines without the header symbol >, loop thru all of them:
  for cross_barcode_seq in $(grep -v "^>" "$cross_barcodes_fasta" ); do
      
     #get line in barcodes fasta with cross sequence
    cross_seq_line=$(sed -n "/${cross_barcode_seq}/=" ${cross_barcodes_fasta});      
    echo "cross_seq_line = ${cross_seq_line}";  
    #get line in barcodes fasta with fasta name
    cross_header_line=$((cross_seq_line -1)); 
    #search reads FASTQ  for the above barcode seqs, -BA get lines above and below match; spit out subset of reads containing the sequence:
    grep -B 1 -A 2 --no-group-separator "${cross_barcode_seq}" reads_with_plaque_barcode.fastq > reads_with_cross_plaque_barcode.fastq;
    echo "grepped the plaque sorted reads fastq on barcode sequence ${cross_barcode_seq} and output matching reads to reads_with_cross_plaque_cross_barcode.fastq";
    num_reads_with_cross_plaque_barcode=$( wc -l reads_with_cross_plaque_barcode.fastq | awk '{print $1/4}')
    #num_reads_with_cross_plaque_barcode=$(grep -c "${cross_barcode_seq}" ${reads_fastq});   
    echo "num reads_with_cross_plaque_barcode = ${num_reads_with_cross_plaque_barcode}";
    echo "total_reads = ${total_reads}";   
    percent_reads_with_cross_plaque_barcode=$(echo "scale=5 ; $num_reads_with_cross_plaque_barcode/$total_reads"| bc);  
    #get header name
    cross_header=$(sed -n "${cross_header_line}p" ${cross_barcodes_fasta});    
    echo "${percent_reads_with_cross_plaque_barcode} fraction of the reads contain an exact match to the barcode for ${cross_header}";   
    #append fraction to array
    barcode_matches_by_plaque+=("${percent_reads_with_cross_plaque_barcode}");   
    echo "{$barcode_matches_by_cross_plaque[@]}";   

    
  
    #create variables for the barcoode sequence + the GSP fwd primer
    plaque_S_primer="CTTTCGTACAACCGAGTAGG";
    plaque_M_primer="CGCTACGGCGGTATTGTC";
    plaque_L_primer="TCGATGTTCAACTACTACGC";
    #search reads files for the above combos
    num_reads_with_S_primer=$(grep -c "${plaque_S_primer}" reads_with_cross_plaque_barcode.fastq);
    num_reads_with_M_primer=$(grep -c "${plaque_M_primer}" reads_with_cross_plaque_barcode.fastq);
    num_reads_with_L_primer=$(grep -c "${plaque_L_primer}" reads_with_cross_plaque_barcode.fastq);
    #print them out
    echo "reads_with_S_primer = ${num_reads_with_S_primer}";
    echo "reads_with_M_primer = ${num_reads_with_M_primer}";
    echo "reads_with_L_primer = ${num_reads_with_L_primer}";
    
    #calculate % of total reads matvhed to the currect barcode, split by segment too
    percent_S_reads_with_cross_plaque_barcode=$(echo "scale=5 ; $num_reads_with_S_primer/$total_reads"| bc);
    percent_M_reads_with_cross_plaque_barcode=$(echo "scale=5 ; $num_reads_with_M_primer/$total_reads"| bc);
    percent_L_reads_with_cross_plaque_barcode=$(echo "scale=5 ; $num_reads_with_L_primer/$total_reads"| bc);
    #print them out
    echo "${percent_S_reads_with_cross_plaque_barcode} of the reads contain an exact S segment match to the barcode for corss ${cross_header} and plaque ${plaque_header}";
    echo "${percent_M_reads_with_cross_plaque_barcode} of the reads contain an exact M segment match to the barcode forcross ${cross_header} and plaque ${plaque_header}";
    echo "${percent_L_reads_with_cross_plaque_barcode} of the reads contain an exact L segment match to the barcode for cross ${cross_header} and plaque ${plaque_header}";
    #append the fraction to an array
    barcode_S_matches_by_cross_plaque+=("${num_reads_with_S_primer}");
    barcode_M_matches_by_cross_plaque+=("${num_reads_with_M_primer}");
    barcode_L_matches_by_cross_plaque+=("${num_reads_with_L_primer}");
    #print the array
    echo "barcode matches by cross and plaque found by grep:"
    echo "${barcode_S_matches_by_cross_plaque[@]}";
    echo "${barcode_M_matches_by_cross_plaque[@]}";
    echo "${barcode_L_matches_by_cross_plaque[@]}";
    
    #get num reads for current barcode in fastq output by dorado demux:
    plaque_num="${plaque_header#*plaque}";
    cross_num="${cross_header#*cross}";

    
    dorado_S_read_count=$( wc -l demux_plaque_S/demux_cross_plaque_S_barcode${plaque_num}/*_cross_S_barcode${cross_num}.fastq | awk '{print $1/4}');
    dorado_M_read_count=$( wc -l demux_plaque_M/demux_cross_plaque_M_barcode${plaque_num}/*_cross_M_barcode${cross_num}.fastq | awk '{print $1/4}');
    dorado_L_read_count=$( wc -l demux_plaque_L/demux_cross_plaque_L_barcode${plaque_num}/*_cross_L_barcode${cross_num}.fastq | awk '{print $1/4}');

    #set to 0 if files dont exist
    if ! test -f demux_plaque_S/demux_cross_plaque_S_barcode${plaque_num}/*_cross_S_barcode${cross_num}.fastq; then
      dorado_S_read_count="0"
    fi;
    if ! test -f demux_plaque_M/demux_cross_plaque_M_barcode${plaque_num}/*_cross_M_barcode${cross_num}.fastq; then
      dorado_M_read_count="0"
    fi;
    if ! test -f demux_plaque_L/demux_cross_plaque_L_barcode${plaque_num}/*_cross_L_barcode${cross_num}.fastq; then
      dorado_L_read_count="0"
    fi;

    
    dorado_cross_plaque_barcode_S_matches_by_plaque+=("${dorado_S_read_count}");
    dorado_cross_plaque_barcode_M_matches_by_plaque+=("${dorado_M_read_count}");
    dorado_cross_plaque_barcode_L_matches_by_plaque+=("${dorado_L_read_count}");
    
    echo "barcode matches by plaque found by dorado:"
    echo "${dorado_cross_plaque_barcode_S_matches_by_plaque[@]}";
    echo "${dorado_cross_plaque_barcode_M_matches_by_plaque[@]}";
    echo "${dorado_cross_plaque_barcode_L_matches_by_plaque[@]}";
    
    echo
    echo "Dorado found ${dorado_S_read_count} reads with the plaque ${plaque_header} and cross ${cross_header} barcode and the S primer | grep found ${num_reads_with_S_primer} reads with the S primer among ${num_reads_with_cross_plaque_barcode} reads with the plaque ${plaque_header} and cross ${cross_header} barcode."
    echo "Dorado found ${dorado_M_read_count} reads with the plaque ${plaque_header} and cross ${cross_header} barcode and the M primer | grep found ${num_reads_with_M_primer} reads with the M primer among ${num_reads_with_cross_plaque_barcode} reads with the plaque ${plaque_header} and cross ${cross_header} barcode."
    echo "Dorado found ${dorado_L_read_count} reads with the plaque ${plaque_header} and cross ${cross_header} barcode and the L primer | grep found ${num_reads_with_L_primer} reads with the L primer among ${num_reads_with_cross_plaque_barcode} reads with the plaque ${plaque_header} and cross ${cross_header} barcode."
   
    echo "${plaque_num},${cross_num},${num_reads_with_plaque_barcode},${dorado_S_read_count},${num_reads_with_S_primer},${dorado_M_read_count},${num_reads_with_M_primer},${dorado_L_read_count},${num_reads_with_L_primer}" >> dorado_vs_grep_barcode_matches.csv;
    echo;   
    echo; 
    done;
  done |& tee update_cross_plaque_barcoding_loop_log.out


  
  done 

      
    done; 
  done; 
done



test
tes
tes


file="dorado_barcode_arrs/plaque_barcodes.fasta"
grep "^>" "$FASTA_FILE" | while IFS= read -r header; do
    # Remove the '>' character from the header
    sequence_id="${header#>"}"

    # You can now use $sequence_id to process each sequence
    echo "Processing sequence: $sequence_id"

done < "${$FASTA_FILE}"

    # Example: Extract the sequence data for the current ID (more complex for multi-line sequences)
    # This example assumes single-line sequences or you need more advanced parsing for multi-line
    # For robust parsing of multi-line sequences, consider tools like 'seqtk subseq' or 'samtools faidx'
    # or using a scripting language like Python with Biopython.
    # grep -A 1 "$header" "$FASTA_FILE" | tail -n 1
done

