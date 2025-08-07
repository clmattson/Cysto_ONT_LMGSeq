

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

file="/path/to/file"
while read line; do
  echo "${line}"
done < "${file}

 
