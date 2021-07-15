#!/bin/bash

### Config
DIS="HCM"      ### Disease (DCM/HCM/ARVC)

### TOOLS
OVERLAP=/hpc/local/CentOS7/dhl_ec/software/overlap.pl

### Directories
MAIN_DIR=/hpc/dhl_ec/aalasiri/Mendelian_CM_MRI_Exome_200K
DIR=${MAIN_DIR}/${DIS}
UKB_200K_WES=/hpc/ukbiobank/WES_200K_2020


### STEP 1 ###

## In ClinVar search ""Arrhythmogenic right ventricular cardiomyopathy"[dis] OR "Arrhythmogenic right ventricular dysplasia"[dis] OR "Arrhythmogenic cardiomyopathy"[dis]" ## put all alternative name of the disease
## Click on pathogenic and Likely pathogenic
## Download file with these features; Format:Tabular (text), Sort by:Location
## Upload file to HPC and rename it with DCM/ARVC/HCM_P_and_LP_select_positionID

### Extract important column
#echo "GRCh38Position Gene dbSNP ID Canonical SPDI Name Review_status" > ${DIR}/${DIS}_P_and_LP_select_positionID
#cat ${DIR}/${DIS}_P_and_LP_v2.txt | sed 's/ /_/g6' | awk -F"\t" '{print $10,$11,$2,$14,$15,$1,$6}' | sed 's/ /:/' >> ${DIR}/${DIS}_P_and_LP_select_positionID

## Extract snps loci  from ClinVar
echo "ID GRCh38Position" > ${DIR}/${DIS}_P_and_LP_positionID
cat ${DIR}/${DIS}_P_and_LP_v3.txt | sed 's/ /_/g6' | awk -F"\t" '{print $15,$10,$11,$15}' | cut -d ' ' -f1,2,3 | sed 's/:/ /1; s/ /:/3; s/_-_/ /g' | awk '{print $2,$3,$2}' | sed 's/:/ /4' | awk '{print $1,$2":"$4}' | tail -n +2 >> ${DIR}/${DIS}_P_and_LP_positionID

## SNP position ID for UKB WES 200K (NO NEED to do it if the file already created)
if [ ! -e ${MAIN_DIR}/ukb_SNP_ID_WES_chrALL.txt ]; then
  for CHR in $(seq 1 22); do
    awk '{print $2, $1":"$4":"$6":"$5}' ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim >> ${MAIN_DIR}/ukb_SNP_ID_WES_chrALL.txt
  done
fi

## Overlap with WES
${OVERLAP} ${DIR}/${DIS}_P_and_LP_positionID 2 ${MAIN_DIR}/ukb_SNP_ID_WES_chrALL.txt 2 | cut -d' ' -f1 >  ${DIR}/overlap_path_SNPs_WES_SNPs.txt
sed 's/:/ /g5' ${MAIN_DIR}/ukb_SNP_ID_WES_chrALL.txt > ${DIR}/ukb_SNP_ID_WES_chrALL.txt_position # replace ":" between position and ref to space
sed 's/:/ /g4' ${DIR}/${DIS}_P_and_LP_positionID > ${DIR}/${DIS}_P_and_LP_positionID_position
${OVERLAP} ${DIR}/${DIS}_P_and_LP_positionID_position 2 ${DIR}/ukb_SNP_ID_WES_chrALL.txt_position 2 > ${DIR}/overlap_path_SNPs_WES_SNPs.txt_position # get overlap by position

### In screen 1
##/hpc/local/CentOS7/dhl_ec/software/overlap.pl overlap_path_SNPs_WES_SNPs.txt 1 overlap_path_SNPs_WES_SNPs.txt_position 1 -v | awk '$1 ~ /I/ {print $0}' | less # get non-overlap Insersion
##/hpc/local/CentOS7/dhl_ec/software/overlap.pl overlap_path_SNPs_WES_SNPs.txt 1 overlap_path_SNPs_WES_SNPs.txt_position 1 -v | awk '$1 ~ /D/ {print $0}' | less  # get non-overlap Deletion
### in screen 2
##less ${DIS}_P_and_LP_positionID ## search in 'less' usin "/" then paste each position from screen 1 "column2", then check if the indel is the same
## compare indel in overlap to the one in ${DIS}_P_and_LP_positionID file, then write them in excel sheet
## emacs overlap_path_SNPs_WES_SNPs.txt # add matched indel snps

### FOR ARVC, combine all SNPs in these two files:
## cat overlapped_ukb_ARVC_all_in_build38_snps >> overlap_path_SNPs_WES_SNPs.txt

### After completing this part, run "run_clinvarSNPs_WES_UKB_script2.sh"
