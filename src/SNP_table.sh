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

echo "Add other important information"
sed 's/:/\t/2' data/temp/CM_incl_SNPs.txt > data/temp/CM_incl
sed -i 's/ /_/g' data/raw/${DIS}_clinvar_result_LP.txt
sed -i 's/ /_/g' data/raw/${DIS}_VKGL.txt
cut -f 10,11,15 data/raw/${DIS}_clinvar_result_LP.txt | sed 's/:/\t/2' | awk '{print $1 ":" $2 ":" $4, $3 ":" $4}' > data/temp/${DIS}_SNP
echo "SNP Canonical_SPDI" > data/temp/${DIS}_SNP_info
bin/overlap.pl data/temp/overlap_SNPs_${DIS}.txt 1 data/temp/${DIS}_SNP 1 >> data/temp/${DIS}_SNP_info

echo "SNP GRCh37 Gene Canonical_SPDI Origin Pathogenicity" | sed 's/ /\t/g' > data/temp/${DIS}_all_SNP
bin/merge_tables.pl --file1 data/raw/${DIS}_clinvar_result_LP.txt --file2 data/temp/${DIS}_SNP_info --index Canonical_SPDI | awk 'NR > 1 {print $1, $10 ":" $11, $4, $2, "ClinVar", $7}' | sed 's/ /\t/g' >> data/temp/${DIS}_all_SNP

cut -f 2-4,9,10,14,16 data/raw/${DIS}_VKGL.txt | awk '{print $2 ":" $4 ":" $5, $1, $3, $6, "VKGL", $7}' > data/temp/${DIS}_VKGL
bin/overlap.pl data/temp/overlap_SNPs_${DIS}.txt 1 data/temp/${DIS}_VKGL 1 > data/temp/${DIS}_VKGL_info
bin/overlap.pl data/temp/${DIS}_all_SNP 1 data/temp/${DIS}_VKGL_info 1 -v | sed 's/ /\t/g' | tail -n +2 >> data/temp/${DIS}_all_SNP
# Manually add the indels to this list using the following commands:
# grep loc data/raw/${DIS}_clinvar_result_LP.txt | awk '{print SNP_38, $8 ":" $9, $2, $15, "ClinVar", $5}' | sed 's/ /\t/g' >> data/temp/${DIS}_all_SNP
# grep loc data/raw/${DIS}_VKGL.txt | awk '{print SNP_38, $2, $4, $14, "VKGL", $16}' | sed 's/ /\t/g' >> data/temp/${DIS}_all_SNP

# Add the Number of individuals per CM to full SNP list
bin/merge_tables.pl --file1 data/temp/overlap_SNPs_ACM.txt --file2 data/temp/CM_incl_SNPs.txt --index SNP > data/temp/tempACM
bin/merge_tables.pl --file1 data/temp/overlap_SNPs_DCM.txt --file2 data/temp/tempACM --index SNP > data/temp/tempDCM
bin/merge_tables.pl --file1 data/temp/overlap_SNPs_HCM.txt --file2 data/temp/tempDCM --index SNP | sed 's/ /\t/g' > data/temp/tempHCM

# Add the extra info to full SNP list
cat data/temp/ACM_all_SNP data/temp/DCM_all_SNP data/temp/HCM_all_SNP | sort -ur > data/temp/All_SNP
bin/merge_tables.pl --file1 data/temp/All_SNP --file2 data/temp/tempHCM --index SNP > data/temp/SNP_table

echo "Add MAF to table"
echo "SNP MAF" > data/temp/SNPs_MAF.frq
cat data/temp/ACM_WES_MRI_UKB_chrALL.allele_frq_filtered data/temp/DCM_WES_MRI_UKB_chrALL.allele_frq_filtered data/temp/HCM_WES_MRI_UKB_chrALL.allele_frq_filtered | awk '{print $2, $5}' | sort -u >> data/temp/SNPs_MAF.frq

bin/merge_tables.pl --file1 data/temp/SNPs_MAF.frq --file2 data/temp/SNP_table --index SNP | sed 's/ /\t/g' > data/processed/overlap_SNPs_CMs.tsv


rm data/temp/CM_incl* data/temp/overlap_SNPs_ACM.txt data/temp/overlap_SNPs_DCM.txt data/temp/overlap_SNPs_HCM.txt data/temp/temp* data/temp/*_VKGL* data/temp/*_SNP_* data/temp/*_SNP data/temp/SNP_table data/temp/SNPs_MAF.frq
