

#run Dorado basecaller on Rosalind:
~/dorado/bin/dorado basecaller -v -r --emit_fastq --disable-read-splitting --no-trim \
-o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall 
sup \
/mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/pod5

#run Dorado demux with custom barcodes.TESTING
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --no-trim --kit-name "custom_barcode_test" \
 --barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/arrangement_test.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/custom_barcode_test.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_test  \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/fastq_pass
 
 
 # DEMULTIPLEX PLAQUES - S
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --no-trim --kit-name "plaque_S" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_plaque_S.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/plaque_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_plaque_S  \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/fastq_pass

 #unclassified : 17333624
#[2025-08-06 16:40:28.445] [debug] Classified rate 14.301747%


 # DEMULTIPLEX PLAQUES - M
 ~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --no-trim --kit-name "plaque_M" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_plaque_M.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/plaque_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_plaque_M  \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/fastq_pass

#[2025-08-06 16:55:20.008] [debug] unclassified : 17511568
#[2025-08-06 16:55:20.008] [debug] Classified rate 13.421983%


 # DEMULTIPLEX PLAQUES - L
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --no-trim --kit-name "plaque_L" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_plaque_L.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/plaque_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_plaque_L  \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/fastq_pass

[2025-08-06 17:19:51.234] [debug] unclassified : 17975868
[2025-08-06 17:19:51.234] [debug] Classified rate 11.126459%


#okay % classified looking sad for the plaque primers. lets get some code that searches for the barcodes and identifies the % of reads containing an exact match

#run the barcode search code to check



#### CROSS Demultiplexing

 # DEMULTIPLEX CROSSES - S
#run in the demux folder
 
for fastq_file in demux_plaque_S/*.fastq; do
plaque_reads="${fastq_file##*/}"; 
plaque_reads="${plaque_reads#*_}"; 
plaque_reads="${plaque_reads%%.*}"; 
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --no-trim --kit-name "cross_S" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_cross_S.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/cross_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_plaque_S/demux_cross_${plaque_reads} \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/${fastq_file};
 done


## # DEMULTIPLEX CROSSES - M

for fastq_file in demux_plaque_M/*.fastq; do
plaque_reads="${fastq_file##*/}"; 
plaque_reads="${plaque_reads#*_}"; 
plaque_reads="${plaque_reads%%.*}"; 
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --no-trim --kit-name "cross_M" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_cross_M.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/cross_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_plaque_M/demux_cross_${plaque_reads} \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/${fastq_file};
 done |& tee demux_cross_M.out



## # DEMULTIPLEX CROSSES - L

for fastq_file in demux_plaque_L/*.fastq; do
plaque_reads="${fastq_file##*/}"; 
plaque_reads="${plaque_reads#*_}"; 
plaque_reads="${plaque_reads%%.*}"; 
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --no-trim --kit-name "cross_L" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_cross_L.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/cross_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_plaque_L/demux_cross_${plaque_reads} \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/${fastq_file};
 done |& tee demux_cross_L.out





#okay, have the files demuxed with decent read counts for each cross/plaque combo

lets concatenate barcode numbers to get unique identifier for each reads file, so can reuse scripts from fall

for fastq_file in demux_plaque_L/demux_cross_plaque_L_barcode*/*cross_L_barcode*.fastq; 
do plaque_barcode_num="${fastq_file#*barcode}"; 
plaque_barcode_num="${plaque_barcode_num%%/*}";
cross_barcode_num="${fastq_file##*barcode}"; 
cross_barcode_num="${cross_barcode_num%%.*}";
segment="${fastq_file%%/demux_cross*}";
segment="${segment##*_}";
echo "plaque${plaque_barcode_num}"; 
echo "cross${cross_barcode_num}"; 
echo "segment: ${segment}";
echo "cp ${fastq_file} to new directory and file named /demux_all/barcode${cross_barcode_num}${plaque_barcode_num}/barcode${cross_barcode_num}${plaque_barcode_num}_${segment}_reads.fastq"
mkdir demux_all/barcode${cross_barcode_num}${plaque_barcode_num};
cp ${fastq_file} demux_all/barcode${cross_barcode_num}${plaque_barcode_num}/barcode${cross_barcode_num}${plaque_barcode_num}_${segment}_reads.fastq
done

#this one:
for fastq_file in demux_plaque_*/demux_cross_plaque_*_barcode*/*cross_*_barcode*.fastq; 
do plaque_barcode_num="${fastq_file#*barcode}"; 
plaque_barcode_num="${plaque_barcode_num%%/*}";
cross_barcode_num="${fastq_file##*barcode}"; 
cross_barcode_num="${cross_barcode_num%%.*}";
segment="${fastq_file%%/demux_cross*}";
segment="${segment##*_}";
echo "plaque${plaque_barcode_num}"; 
echo "cross${cross_barcode_num}"; 
echo "segment: ${segment}";
echo "cp ${fastq_file} to new directory and file named /demux_all/barcode${cross_barcode_num}${plaque_barcode_num}/barcode${cross_barcode_num}${plaque_barcode_num}_${segment}_reads.fastq"
mkdir demux_all/barcode${cross_barcode_num}${plaque_barcode_num};
cp ${fastq_file} demux_all/barcode${cross_barcode_num}${plaque_barcode_num}/barcode${cross_barcode_num}${plaque_barcode_num}_${segment}_reads.fastq
done




move them to franklin







 #### CROSS Demultiplexing WITH TRIM ENABLED
 
 # DEMULTIPLEX PLAQUES - S
#run in the demux folder


 # DEMULTIPLEX PLAQUES - S
 
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --kit-name "plaque_S" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_plaque_S.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/plaque_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_trimmed_plaque_S  \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/fastq_pass

# M 
 ~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --kit-name "plaque_M" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_plaque_M.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/plaque_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_trimmed_plaque_M  \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/fastq_pass

 ~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --kit-name "plaque_L" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_plaque_L.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/plaque_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_trimmed_plaque_L  \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/fastq_pass


 
#### CROSS Demultiplexing

 # DEMULTIPLEX CROSSES - S
#run in the demux folder
 
for fastq_file in demux_trimmed_plaque_S/*.fastq; do
plaque_reads="${fastq_file##*/}"; 
plaque_reads="${plaque_reads#*_}"; 
plaque_reads="${plaque_reads%%.*}"; 
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --kit-name "cross_S" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_cross_S.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/cross_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_trimmed_plaque_S/demux_trimmed_cross_${plaque_reads} \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/${fastq_file};
 done


## # DEMULTIPLEX CROSSES - M

for fastq_file in demux_trimmed_plaque_M/*.fastq; do
plaque_reads="${fastq_file##*/}"; 
plaque_reads="${plaque_reads#*_}"; 
plaque_reads="${plaque_reads%%.*}"; 
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --kit-name "cross_M" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_cross_M.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/cross_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_trimmed_plaque_M/demux_trimmed_cross_${plaque_reads} \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/${fastq_file};
 done |& tee demux_trimmed_cross_M.out



## # DEMULTIPLEX CROSSES - L

for fastq_file in demux_trimmed_plaque_L/*.fastq; do
plaque_reads="${fastq_file##*/}"; 
plaque_reads="${plaque_reads#*_}"; 
plaque_reads="${plaque_reads%%.*}"; 
~/dorado/bin/dorado demux -v --emit-fastq --emit-summary --kit-name "cross_L" \
--barcode-arrangement /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/arr_cross_L.txt \
 --barcode-sequences /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/dorado_barcode_arrs/cross_barcodes.fasta \
 -o /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/demux_trimmed_plaque_L/demux_trimmed_cross_${plaque_reads} \
 /mnt/data0/MinION_reads/Cysto_LMGSeq_barcode_test_07302025/Cysto_LMGSeq_barcode_test_07302025/20250730_1606_MN23913_FBC73506_d5fa85e0/sup_basecall/${fastq_file};
 done |& tee demux_trimmed_cross_L.out


#trim primers? ... no, unnecessary since part of sequence


#redo for trimmed data
for fastq_file in demux_trimmed_plaque_*/demux_trimmed_cross_plaque_*_barcode*/*cross_*_barcode*.fastq; 
do plaque_barcode_num="${fastq_file#*barcode}"; 
plaque_barcode_num="${plaque_barcode_num%%/*}";
cross_barcode_num="${fastq_file##*barcode}"; 
cross_barcode_num="${cross_barcode_num%%.*}";
segment="${fastq_file%%/demux_trimmed_cross*}";
segment="${segment##*_}";
echo "plaque${plaque_barcode_num}"; 
echo "cross${cross_barcode_num}"; 
echo "segment: ${segment}";
echo "cp ${fastq_file} to new directory and file named /demux_all_trimmed/barcode${cross_barcode_num}${plaque_barcode_num}/barcode${cross_barcode_num}${plaque_barcode_num}_${segment}_reads.fastq"
mkdir demux_all_trimmed/barcode${cross_barcode_num}${plaque_barcode_num};
cp ${fastq_file} demux_all_trimmed/barcode${cross_barcode_num}${plaque_barcode_num}/barcode${cross_barcode_num}${plaque_barcode_num}_${segment}_reads.fastq
done









#running scripts
bash make_databases.sh -d /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all -e /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/sample_list2.csv -c /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/cross_list.txt -s /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/refs/ref_phi6_S_04.fasta -m /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/refs/ref_phi6_M_01.fasta -l /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/refs/ref_phi6_L_89.fasta


bash make_databases.sh -d /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all -e /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/sample_list2.csv -c /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/cross_list.txt -s /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/refs/ref_phi6_S_04.fasta -m /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/refs/ref_phi6_M_01.fasta -l /group/sldmunozgrp/cysto_LMGSeq08-25/demux_all/refs/ref_phi6_L_89.fasta




#cutadapt
bash make_databases.sh -d /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt -e /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/sample_list_plate01.csv -c /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/cross_list.txt -s /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_S_04.fasta -m /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_M_01.fasta -l /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_L_89.fasta


bash genotyping_new.sh -d /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt -e /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/sample_list_plate01.csv -c /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/cross_list.txt -s /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_S_04.fasta -m /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_M_01.fasta -l /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/old_ref_phi6_L_89.fasta


bash genotyping_new.sh -d /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt -e /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/sample_list2.csv -c /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/cross_list.txt -s /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_S_04.fasta -m /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_M_01.fasta -l /group/sldmunozgrp/cysto_LMGSeq08-25/cutadapt/refs/ref_phi6_L_89.fasta
