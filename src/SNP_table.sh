#!/bin/bash

echo "Make table of all SNPs with their respective numbers per CM"
echo "SNP" > data/temp/CM_incl_SNPs.txt
cat data/temp/ACM_overlap_LP_WES_SNPs_updated.txt data/temp/DCM_overlap_LP_WES_SNPs_updated.txt data/temp/HCM_overlap_LP_WES_SNPs_updated.txt | sort -u >> data/temp/CM_incl_SNPs.txt

echo "SNP N_${DIS}" > data/temp/overlap_SNPs_${DIS}.txt
while IFS= read -r line; do

  echo "Searching for the following SNP now: ${line}"
  COL=`head -1 data/temp/${DIS}_temp/${DIS}_SNPs_MRI.txt | bin/transpose_perl.pl | grep -n ${line} | sed 's/:/\t/' | cut -f1`
  echo -e "${line} $(awk -v col="${COL}" '{print $col}' data/temp/${DIS}_temp/${DIS}_SNPs_MRI.txt | grep -c 1)" >> data/temp/overlap_SNPs_${DIS}.txt

done < "data/temp/${DIS}_overlap_LP_WES_SNPs_updated.txt"

bin/merge_tables.pl --file1 data/temp/overlap_SNPs_ACM.txt --file2 data/temp/CM_incl_SNPs.txt --index SNP > data/temp/temp1
bin/merge_tables.pl --file1 data/temp/overlap_SNPs_DCM.txt --file2 data/temp/temp1 --index SNP > data/temp/temp2
bin/merge_tables.pl --file1 data/temp/overlap_SNPs_HCM.txt --file2 data/temp/temp2 --index SNP | sed 's/ /\t/g' > data/temp/temp3


echo "Add MAF to table"
echo "SNP MAF" > data/temp/SNPs_MAF.frq
cat data/temp/ACM_WES_MRI_UKB_chrALL.allele_frq_filtered data/temp/DCM_WES_MRI_UKB_chrALL.allele_frq_filtered data/temp/HCM_WES_MRI_UKB_chrALL.allele_frq_filtered | awk '{print $2, $5}' | sort -u >> data/temp/SNPs_MAF.frq

bin/merge_tables.pl --file1 data/temp/SNPs_MAF.frq --file2 data/temp/temp3 --index SNP | sed 's/ /\t/g' > data/temp/overlap_SNPs_CMs.tsv

rm data/temp/CM_incl_SNPs.txt data/temp/overlap_SNPs_ACM.txt data/temp/overlap_SNPs_DCM.txt data/temp/overlap_SNPs_HCM.txt data/temp/temp*
