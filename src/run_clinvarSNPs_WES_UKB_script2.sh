#!/bin/bash
### TOOLS
PLINK=/hpc/local/CentOS7/dhl_ec/software/plink_v1.9
OVERLAP=/hpc/local/CentOS7/dhl_ec/software/overlap.pl
TRNASPOSE=/hpc/dhl_ec/aalasiri/ANGPTL3/transpose_perl.pl

DIS="HCM"
MRI="LV"
MAIN_DIR=/hpc/dhl_ec/aalasiri/Mendelian_CM_MRI_Exome_200K/
DIR=/hpc/dhl_ec/aalasiri/Mendelian_CM_MRI_Exome_200K/${DIS}
UKB_200K_WES=/hpc/ukbiobank/WES_200K_2020
UKB_200K_WES_FAM=/hpc/dhl_ec/data/ukbiobank/WES_200K_2020

### STEP 2 ###
### overlap samples betwen 200K WES and LV
## get IID of LV/RV MRI after QC
#tail -n +2 /hpc/dhl_ec/aalasiri/CMR_metrics/PheWAS_MRI/ukb_MRI_LV_clean.txt | awk '{print $1,$1}' > ${DIR}/ukb_MRI_LV_IID_clean.txt
tail -n +2 /hpc/dhl_ec/aalasiri/CMR_metrics/PheWAS_MRI/ukb_MRI_${MRI}_clean.txt | awk '{print $1,$1}' > ${DIR}/ukb_MRI_${MRI}_IID_clean.txt
## Allele freq. using Plink
rm ${DIR}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq
for CHR in $(seq 1 22); do
  echo $CHR
  $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
         --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
         --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
         --keep ${DIR}/ukb_MRI_${MRI}_IID_clean.txt \
         --extract ${DIR}/overlap_path_SNPs_WES_SNPs.txt \
         --freq --out ${DIR}/${DIS}_WES_MRI_UKB_chr${CHR}
  tail -n +2 ${DIR}/${DIS}_WES_MRI_UKB_chr${CHR}.frq >> ${DIR}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq
  rm ${DIR}/${DIS}_WES_MRI_UKB_chr${CHR}.*
done

awk '$5 != 0 {print $0}' ${DIR}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq > ${DIR}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq_filtered
awk '{print $2}' ${DIR}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq_filtered > ${DIR}/${DIS}_${MRI}_overlap_path_SNPs_WES_SNPs_updated.txt

## Extract filtered SNPs and samples in VCF file For LV sample size
#for CHR in $(seq 1 22); do
#  $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
#         --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
#         --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
#         --keep ${DIR}/ukb_MRI_${MRI}_IID_clean.txt \
#         --extract ${DIR}/${DIS}_${MRI}_overlap_path_SNPs_WES_SNPs_updated.txt \
#         --recode vcf-iid \
#         --keep-allele-order --out ${DIR}/${DIS}_WES_MRI_UKB_chr${CHR}
#done

#HEADER_CHR=$(head -1 ${DIR}/${DIS}_overlap_path_SNPs_WES_SNPs_updated.txt | cut -d':' -f1)
#head -7 ${DIR}/${DIS}_WES_MRI_UKB_chr${HEADER_CHR}.vcf > ${DIR}/${DIS}_WES_MRI_UKB_chrALL.vcf
#for CHR in $(seq 1 22); do
#  tail -n +8 ${DIR}/${DIS}_WES_MRI_UKB_chr${CHR}.vcf >> ${DIR}/${DIS}_WES_MRI_UKB_chrALL.vcf
#  rm ${DIR}/${DIS}_WES_MRI_UKB_chr${CHR}.*
#done
#tail -n +7 ${DIR}/${DIS}_WES_MRI_UKB_chrALL.vcf  | cut -f3,10- | ${TRNASPOSE} > ${DIR}/${DIS}_WES_MRI_UKB_chrALL.transpose

## Get samples ID for each SNP carrier
#SNP_LIST=$(head -1 ${DIR}/${DIS}_WES_MRI_UKB_chrALL.transpose | sed 's/\t/ /g' | cut -d' ' -f2-)
#rm ${DIR}/${DIS}_WES_MRI_UKB_snps_IID.txt
#for SNP in $SNP_LIST; do
#  echo $SNP
#  echo $SNP >> ${DIR}/${DIS}_WES_MRI_UKB_snps_IID.txt
#  awk -v VAR=$SNP 'NR==1{for(i=1; i<=NF; i++) if ($i==VAR) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf $1 "\n"}' ${DIR}/${DIS}_WES_MRI_UKB_chrALL.transpose | awk '$1 != "0/0" {print $2,$1}' | tail -n +2 | cut -d' ' -f1 >> ${DIR}/${DIS}_WES_MRI_UKB_snps_IID.txt
#done

### STEP 3 ###
### For ALL WES samples
rm ${DIR}/WES_MRI_UKB_chrALL.allele_frq
for CHR in $(seq 1 22); do
  $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
         --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
         --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
         --extract ${DIR}/overlap_path_SNPs_WES_SNPs.txt \
         --freq --out ${DIR}/WES_MRI_UKB_chr${CHR}
  tail -n +2 ${DIR}/WES_MRI_UKB_chr${CHR}.frq >> ${DIR}/WES_MRI_UKB_chrALL.allele_frq
  rm ${DIR}/WES_MRI_UKB_chr${CHR}.*
done

## Extract filtered SNPs and samples in VCF file For all 200K WES
#awk '$5 != 0 {print $0}' ${DIR}/WES_MRI_UKB_chrALL.allele_frq > ${DIR}/WES_MRI_UKB_chrALL.allele_frq_filtered
#awk '{print $2}' ${DIR}/WES_MRI_UKB_chrALL.allele_frq_filtered > ${DIR}/overlap_path_SNPs_WES_SNPs_updated.txt
#for CHR in $(seq 1 22); do
#  $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
#         --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
#         --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
#         --extract ${DIR}/overlap_path_SNPs_WES_SNPs_updated.txt \
#         --recode vcf-iid \
#         --keep-allele-order --out ${DIR}/WES_MRI_UKB_chr${CHR}
#done
#HEADER_CHR_2=$(head -1 ${DIR}/overlap_path_SNPs_WES_SNPs_updated.txt | cut -d':' -f1)
#head -7 ${DIR}/WES_MRI_UKB_chr${HEADER_CHR_2}.vcf > ${DIR}/WES_MRI_UKB_chrALL.vcf
#for CHR in $(seq 1 22); do
#  tail -n +8 ${DIR}/WES_MRI_UKB_chr${CHR}.vcf >> ${DIR}/WES_MRI_UKB_chrALL.vcf
#  rm ${DIR}/WES_MRI_UKB_chr${CHR}.*
#done
#tail -n +7 ${DIR}/WES_MRI_UKB_chrALL.vcf  | cut -f3,10- | ${TRNASPOSE} > ${DIR}/WES_MRI_UKB_chrALL.transpose
#SNP_LIST_2=$(head -1 ${DIR}/WES_MRI_UKB_chrALL.transpose | sed 's/\t/ /g' | cut -d' ' -f2-)
## Get samples ID for each SNP carrier
#rm ${DIR}/WES_MRI_UKB_snps_IID.txt
#for SNP in $SNP_LIST_2; do
#  echo $SNP
#  echo $SNP >> ${DIR}/WES_MRI_UKB_snps_IID.txt
#  awk -v VAR=$SNP 'NR==1{for(i=1; i<=NF; i++) if ($i==VAR) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf $1 "\n"}' ${DIR}/WES_MRI_UKB_chrALL.transpose | awk '$1 != "0/0" {print $2,$1}' | tail -n +2 | cut -d' ' -f1 >> ${DIR}/WES_MRI_UKB_snps_IID.txt
#done
