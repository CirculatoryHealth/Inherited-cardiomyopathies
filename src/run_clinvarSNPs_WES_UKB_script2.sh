#! /bin/bash
set -e

### Creating display functions
### Setting colouring
NONE='\033[00m'
OPAQUE='\033[2m'
FLASHING='\033[5m'
BOLD='\033[1m'
ITALIC='\033[3m'
UNDERLINE='\033[4m'
STRIKETHROUGH='\033[9m'

RED='\033[01;31m'
GREEN='\033[01;32m'
YELLOW='\033[01;33m'
PURPLE='\033[01;35m'
CYAN='\033[01;36m'
WHITE='\033[01;37m'

function echobold { #'echobold' is the function name
    echo -e "${BOLD}${1}${NONE}" # this is whatever the function needs to execute, note ${1} is the text for echo
}
function echoitalic {
    echo -e "${ITALIC}${1}${NONE}"
}
function echonooption {
    echo -e "${OPAQUE}${RED}${1}${NONE}"
}
function echoerrorflash {
    echo -e "${RED}${BOLD}${FLASHING}${1}${NONE}"
}
function echoerror {
    echo -e "${RED}${1}${NONE}"
}
# errors no option
function echoerrornooption {
    echo -e "${YELLOW}${1}${NONE}"
}
function echoerrorflashnooption {
    echo -e "${YELLOW}${BOLD}${FLASHING}${1}${NONE}"
}
function importantnote {
    echo -e "${CYAN}${1}${NONE}"
}

script_copyright_message() {
	echo ""
	THISYEAR=$(date +'%Y')
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "+                                                                                                       +"
	echo "+ Copyright (c) 2021-${THISYEAR} Abdulrahman Alasiri                                                           +"
	echo "+                                                                                                       +"
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
}

script_arguments_error() {
	echoerror "$1" # Additional message
	echoerror "- Argument #1 is selecting [ARVC], [DCM] or [HCM]."
  echoerror "- Argument #2 is selecting [RV] or [LV]."
	echoerror ""
	echoerror "An example command would be: run_clinvarSNPs_WES_UKB_script2.sh.sh [arg1] [arg2]"
	echoerror "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
 	echo ""
	script_copyright_message
	exit 1
}

echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "+                Extraction of cardiomyopathies mutations from exomes in UK Biobank                     +"
echobold "+                                                                                                       +"
echobold "+                                                                                                       +"
echobold "+ * Written by  : Abdulrahman Alasiri                                                                   +"
echobold "+ * E-mail      : a.i.alasiri@umcutrecht.nl                                                             +"
echobold "+ * Last update : 2021-02-22                                                                            +"
echobold "+ * Version     : 0.1.0                                                                                 +"
echobold "+                                                                                                       +"
echobold "+ * Description : This script will extract [Likely] pathogenic mutations                                +"
echobold "+                 for ARVC, DCM and HCM from UKBB exomes PLINK files.                                   +"
echobold "+                                                                                                       +"
echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's date and time: "$(date)
TODAY=$(date +"%Y%m%d")
echo ""

### START of if-else statement for the number of command-line arguments passed ###
if [[ $# -lt 2 ]]; then
    echoerrorflash "                                     *** Oh no! Computer says no! ***"
    echo ""
    script_arguments_error "You must supply at least [2] argument when running extraction of cardiomyopathy mutations from 200K WES in UKBB!"

else
    ### REQUIRED | GENERALS
    DIS="$1"    # Disease (DCM/HCM/ARVC)
    MRI="$2"    # LV or RV
    ### TOOLS
    PLINK=/hpc/local/CentOS7/dhl_ec/software/plink_v1.9
    OVERLAP=/hpc/local/CentOS7/dhl_ec/software/overlap.pl
    TRNASPOSE=/hpc/dhl_ec/aalasiri/scripts/Inherited_cardiomyopathies/bin/transpose_perl.pl

    ### PROJECT SPECIFIC
    PROJECTDIR=/hpc/dhl_ec/aalasiri/scripts/Inherited_cardiomyopathies # the root directoryto repo
    INPUT=${PROJECTDIR}/data/raw
    OUTPUT=${PROJECTDIR}/data/processed
    TEMP=${PROJECTDIR}/data/temp
    UKB_200K_WES=/hpc/ukbiobank/WES_200K_2020
    UKB_200K_WES_FAM=/hpc/dhl_ec/data/ukbiobank/WES_200K_2020


### STEP 2 ###
### overlap samples betwen 200K WES and LV
## get IID of LV/RV MRI after QC
    tail -n +2 /hpc/dhl_ec/aalasiri/CMR_metrics/PheWAS_MRI/ukb_MRI_${MRI}_clean.txt | awk '{print $1,$1}' > ${TEMP}/ukb_MRI_${MRI}_IID_clean.txt
## Allele freq. using Plink
    rm ${OUTPUT}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq
    for CHR in $(seq 1 22) X Y; do
      echo $CHR
      $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
             --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
             --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
             --keep ${TEMP}/ukb_MRI_${MRI}_IID_clean.txt \
             --extract ${TEMP}/${DIS}_overlap_path_SNPs_WES_SNPs.txt \
             --freq --out ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}
      tail -n +2 ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.frq >> ${OUTPUT}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq
      rm ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.*
    done

    awk '$5 != 0 {print $0}' ${OUTPUT}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq > ${OUTPUT}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq_filtered
    awk '{print $2}' ${OUTPUT}/${DIS}_WES_${MRI}_MRI_UKB_chrALL.allele_frq_filtered > ${OUTPUT}/${DIS}_${MRI}_overlap_path_SNPs_WES_SNPs_updated.txt

## Extract filtered SNPs and samples in VCF file For LV sample size
    for CHR in $(seq 1 22) X Y; do
      $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
             --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
             --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
             --keep ${TEMP}/ukb_MRI_${MRI}_IID_clean.txt \
             --extract ${OUTPUT}/${DIS}_${MRI}_overlap_path_SNPs_WES_SNPs_updated.txt \
             --recode vcf-iid \
             --keep-allele-order --out ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}
    done

    HEADER_CHR=$(head -1 ${OUTPUT}/${DIS}_${MRI}_overlap_path_SNPs_WES_SNPs_updated.txt | cut -d':' -f1)
    head -7 ${TEMP}/${DIS}_WES_MRI_UKB_chr${HEADER_CHR}.vcf > ${OUTPUT}/${DIS}_WES_MRI_UKB_chrALL.vcf
    for CHR in $(seq 1 22); do
      tail -n +8 ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.vcf >> ${OUTPUT}/${DIS}_WES_MRI_UKB_chrALL.vcf
      rm ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.*
    done
    tail -n +7 ${OUTPUT}/${DIS}_WES_MRI_UKB_chrALL.vcf  | cut -f3,10- | ${TRNASPOSE} > ${OUTPUT}/${DIS}_WES_MRI_UKB_chrALL.transpose

    ## Get samples ID for each SNP carrier
    SNP_LIST=$(head -1 ${OUTPUT}/${DIS}_WES_MRI_UKB_chrALL.transpose | sed 's/\t/ /g' | cut -d' ' -f2-)
    rm ${OUTPUT}/${DIS}_WES_MRI_UKB_snps_IID.txt
    for SNP in $SNP_LIST; do
      echo $SNP
      echo $SNP >> ${OUTPUT}/${DIS}_WES_MRI_UKB_snps_IID.txt
      awk -v VAR=$SNP 'NR==1{for(i=1; i<=NF; i++) if ($i==VAR) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf $1 "\n"}' ${OUTPUT}/${DIS}_WES_MRI_UKB_chrALL.transpose | awk '$1 != "0/0" {print $2,$1}' | tail -n +2 | cut -d' ' -f1 >> ${OUTPUT}/${DIS}_WES_MRI_UKB_snps_IID.txt
    done

### STEP 3 ###
### For ALL 200K WES samples
    rm ${OUTPUT}/WES_MRI_UKB_chrALL.allele_frq
    for CHR in $(seq 1 22); do
      $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
             --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
             --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
             --extract ${TEMP}/${DIS}_overlap_path_SNPs_WES_SNPs.txt \
             --freq --out ${TEMP}/WES_MRI_UKB_chr${CHR}
      tail -n +2 ${TEMP}/WES_MRI_UKB_chr${CHR}.frq >> ${OUTPUT}/WES_MRI_UKB_chrALL.allele_frq
      rm ${TEMP}/WES_MRI_UKB_chr${CHR}.*
    done

## Extract filtered SNPs and samples in VCF file For all 200K WES
    awk '$5 != 0 {print $0}' ${OUTPUT}/WES_MRI_UKB_chrALL.allele_frq > ${OUTPUT}/WES_MRI_UKB_chrALL.allele_frq_filtered
    awk '{print $2}' ${OUTPUT}/WES_MRI_UKB_chrALL.allele_frq_filtered > ${OUTPUT}/overlap_path_SNPs_WES_SNPs_updated.txt
    for CHR in $(seq 1 22); do
      $PLINK --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
             --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
             --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
             --extract ${OUTPUT}/${DIS}_${MRI}_overlap_path_SNPs_WES_SNPs_updated.txt \
             --recode vcf-iid \
             --keep-allele-order --out ${TEMP}/WES_MRI_UKB_chr${CHR}
    done
    HEADER_CHR_2=$(head -1 ${OUTPUT}/${DIS}_${MRI}_overlap_path_SNPs_WES_SNPs_updated.txt  | cut -d':' -f1)
    head -7 ${TEMP}/WES_MRI_UKB_chr${HEADER_CHR_2}.vcf > ${OUTPUT}/WES_MRI_UKB_chrALL.vcf
    for CHR in $(seq 1 22); do
      tail -n +8 ${TEMP}/WES_MRI_UKB_chr${CHR}.vcf >> ${OUTPUT}/WES_MRI_UKB_chrALL.vcf
      rm ${TEMP}/WES_MRI_UKB_chr${CHR}.*
    done
    tail -n +7 ${OUTPUT}/WES_MRI_UKB_chrALL.vcf  | cut -f3,10- | ${TRNASPOSE} > ${OUTPUT}/WES_MRI_UKB_chrALL.transpose
    SNP_LIST_2=$(head -1 ${OUTPUT}/WES_MRI_UKB_chrALL.transpose | sed 's/\t/ /g' | cut -d' ' -f2-)
# Get samples ID for each SNP carrier
    rm ${OUTPUT}/WES_MRI_UKB_snps_IID.txt
    for SNP in $SNP_LIST_2; do
      echo $SNP
      echo $SNP >> ${OUTPUT}/WES_MRI_UKB_snps_IID.txt
      awk -v VAR=$SNP 'NR==1{for(i=1; i<=NF; i++) if ($i==VAR) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf $1 "\n"}' ${OUTPUT}/WES_MRI_UKB_chrALL.transpose | awk '$1 != "0/0" {print $2,$1}' | tail -n +2 | cut -d' ' -f1 >> ${OUTPUT}/WES_MRI_UKB_snps_IID.txt
    done
