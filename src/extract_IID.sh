#!/bin/bash
#
#SBATCH --job-name extract_IID                                             		  # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/extract_IID.log              # log file of this script
#SBATCH --mem=10G                                                          		  # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                                   # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                                        # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                            Script to extract individuals carrying (likely) pathogenic CM mutations                           "
echo "                                                    version 1.0 (20210714)                                                    "
echo ""
echo "* Written by      : Marion van Vugt"
echo "* E-mail          : m.vanvugt-2@umcutrecht.nl"
echo "* Last update     : 2021-07-14"
echo "* Version         : extract_IID_1.0"
echo ""
echo "* Description     : This script will extract individuals carrying (likely) pathogenic CM mutations     "
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

echo ""
date "+DATE: %a %d/%m/%Y%nTIME: %T"

echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

echo ""

NONE='\e[00m'
BOLD='\e[1m'
RED='\e[01;31m'


function echobold {                       # 'echobold' is the function name
    echo -e "${BOLD}${1}${NONE}"          # this is whatever the function needs to execute, note ${1} is the text for echo
}

function echoerror {
    echo -e "${RED}${1}${NONE}"
}

script_copyright_message() {
	echo ""
	THISYEAR=$(date +'%Y')
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "+ The MIT License (MIT)                                                                                 +"
	echo "+ Copyright (c) 1979-${THISYEAR} Marion van Vugt                                                        +"
	echo "+                                                                                                       +"
	echo "+ Permission is hereby granted, free of charge, to any person obtaining a copy of this software and     +"
	echo "+ associated documentation files (the \"Software\"), to deal in the Software without restriction,       +"
	echo "+ including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, +"
	echo "+ and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, +"
	echo "+ subject to the following conditions:                                                                  +"
	echo "+                                                                                                       +"
	echo "+ The above copyright notice and this permission notice shall be included in all copies or substantial  +"
	echo "+ portions of the Software.                                                                             +"
	echo "+                                                                                                       +"
	echo "+ THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT   +"
	echo "+ NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                +"
	echo "+ NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES  +"
	echo "+ OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN   +"
	echo "+ CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                            +"
	echo "+                                                                                                       +"
	echo "+ Reference: http://opensource.org.                                                                     +"
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
}

script_arguments_error() {
	echoerror "$1" # ERROR MESSAGE
  echoerror "- Argument #1   --  Disease to perform analysis for, could be 'DCM'"
  echoerror "- Argument #2   --  File name and path of/to the overlapping indels, could be '/hpc/dhl_ec/mvanvugt/UKBB/indels.txt'"
  echoerror "- Argument #3   --  File name and path of/to the log file, could be '/hpc/dhl_ec/mvanvugt/UKBB/log.txt'"
	echoerror ""
	echoerror "An example command would be: extract_IID.sh [arg1: DCM] [arg2: /hpc/dhl_ec/mvanvugt/UKBB/indels.txt] [arg3: /hpc/dhl_ec/mvanvugt/UKBB/log.txt]."
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  	# The wrong arguments are passed, so we'll exit the script now!
  	exit 1
}

if [[ $# -lt 3 ]]; then
  echo "Error, number of arguments found "$#"."
  script_arguments_error "You must supply [3] correct arguments when running this script"

else

  DIS="$1"      ### Disease (DCM/HCM/ACM)
  INDEL="$2"
  LOG="$3"
  DIR="data/temp"
  ROOT=$(pwd)

  if [[ ! -d ${DIR}/temp ]]; then
    mkdir -v ${DIR}/temp
  fi
  TEMP="${DIR}/temp"
  PLINK="/hpc/local/CentOS7/dhl_ec/software/plink_v1.9"
  UKB_200K_WES="/hpc/ukbiobank/WES_200K_2020"
  UKB_200K_WES_FAM='/hpc/dhl_ec/data/ukbiobank/WES_200K_2020'


  echo ""
  echo "Settings for this script:"
  echo "Disease: ____________________ [ ${DIS} ]"
  echo "Directory: __________________ [ ${DIR} ]"
  echo "Indel file: _________________ [ ${INDEL} ]"
  echo "Temporary directory: ________ [ ${TEMP} ]"

  echo ""

  ### STEP 2 ###
  ### overlap samples betwen 200K WES and LV/RV
  ## get IID of LV/RV MRI after QC
  tail -n +2 /hpc/dhl_ec/aalasiri/CMR_metrics/PheWAS_MRI/ukb_MRI_LV_clean.txt | awk '{print $1,$1}' > ${TEMP}/ukb_MRI_LV_IID_clean.txt
  tail -n +2 /hpc/dhl_ec/aalasiri/CMR_metrics/PheWAS_MRI/ukb_MRI_RV_clean.txt | awk '{print $1,$1}' > ${TEMP}/ukb_MRI_RV_IID_clean.txt

  # Adding indels to SNP-lists
  echo -e "$(wc -l ${INDEL} | cut -d" " -f1) indels included" >> ${TEMP}/${DIS}_log
  awk '{print $1}' ${INDEL} >> ${DIR}/${DIS}_overlap_LP_WES_SNPs.txt

  ## Allele freq. using Plink
  rm ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.allele_frq
  for CHR in $(seq 1 22); do
    echo ${CHR}
    ${PLINK} --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
    --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
    --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
    --extract ${DIR}/${DIS}_overlap_LP_WES_SNPs.txt \
    --freq --out ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}
    tail -n +2 ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.frq >> ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.allele_frq
    echo ""
  done

  awk '$5 != 0 {print $0}' ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.allele_frq > ${DIR}/${DIS}_WES_MRI_UKB_chrALL.allele_frq_filtered
  awk '{print $2}' ${DIR}/${DIS}_WES_MRI_UKB_chrALL.allele_frq_filtered > ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt
  echo -e "$(wc -l ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt | cut -d" " -f1) SNPs with carriers in the UK Biobank" >> ${TEMP}/${DIS}_log

  # Remove SNPs that have overlap between HCM/DCM/ACM and are not annotated well
  if [[ -e data/raw/${DIS}_remove_overlap.txt ]]; then
    echo -e "$(wc -l data/raw/${DIS}_remove_overlap.txt | cut -d" " -f1) SNPs removed due to ambiguous annotation" >> ${TEMP}/${DIS}_log
    awk 'NR==FNR{a[$0];next}!($0 in a)' data/raw/${DIS}_remove_overlap.txt ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt > ${TEMP}/${DIS}_temp
    mv ${TEMP}/${DIS}_temp ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt
  fi

  echo -e "$(wc -l ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt | cut -d" " -f1) SNPs to extract UKB-individuals for" >> ${TEMP}/${DIS}_log
  echo "Making vcf-files with only desired snps"
  for CHR in $(seq 1 22); do
    echo ${CHR}
    ${PLINK} --bed ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bed \
    --bim ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim \
    --fam ${UKB_200K_WES_FAM}/UKBexomeOQFE_200K_chr${CHR}_v1.fam \
    --extract ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt \
    --recode vcf-iid \
    --keep-allele-order --out ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}
    echo ""
  done

  echo "Making one vcf-file"
  HEADER_CHR_2=$(head -1 ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt | cut -d':' -f1)
  head -7 ${TEMP}/${DIS}_WES_MRI_UKB_chr${HEADER_CHR_2}.vcf > ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.vcf
  for CHR in $(seq 1 22); do
    tail -n +8 ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.vcf >> ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.vcf
  done

  echo "Transposing vcf-file so header are SNPs and rows are individuals"
  tail -n +7 ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.vcf  | cut -f3,10- | bin/transpose_perl.pl > ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.transpose

  # Get number of lines
  NUM=`awk 'NR == 1 {print NF}' ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.transpose`
  echo ""
  echo "Making file with individuals per SNP"
  for SNP in $(seq 2 ${NUM}); do
    cut -f 1,${SNP} ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.transpose | awk '$2 != "0/0" && $2 != "./." {print}' > ${TEMP}/${DIS}_SNP${SNP}_ID.txt
  done

  # Make a list of all individuals
  echo "ID" > ${TEMP}/${DIS}_LP_allIDs.txt
  cat ${TEMP}/${DIS}_SNP*txt | grep -v : | awk '{print $1}' | sort -u >> ${TEMP}/${DIS}_LP_allIDs.txt

  echo "Merging all seperate SNP-files"
  for SNP in $(seq 2 ${NUM}); do
    if [[ ${SNP} == 2 ]]; then
      bin/merge_tables.pl --file1 ${TEMP}/${DIS}_SNP${SNP}_ID.txt --file2 ${TEMP}/${DIS}_LP_allIDs.txt --index ID > ${TEMP}/${DIS}_merge${SNP}
    else
      MIN=$((${SNP}-1))
      bin/merge_tables.pl --file1 ${TEMP}/${DIS}_SNP${SNP}_ID.txt --file2 ${TEMP}/${DIS}_merge${MIN} --index ID > ${TEMP}/${DIS}_merge${SNP}
    fi
  done

  echo "Making list of genes per SNP"
  head -1 ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.transpose | bin/transpose_perl.pl > ${TEMP}/${DIS}_snps_temp
  echo "ID Gene" > ${TEMP}/${DIS}_sng_temp
  tail -n +2 ${TEMP}/${DIS}_LP_positionID | awk '{print $2, $4}' >> ${TEMP}/${DIS}_sng_temp
  cat ${INDEL} >> ${TEMP}/${DIS}_sng_temp
  bin/merge_tables.pl --file1 ${TEMP}/${DIS}_sng_temp --file2 ${TEMP}/${DIS}_snps_temp --index ID | awk '{print $2}' > ${TEMP}/${DIS}_genes

  sed 's/0\///g' ${TEMP}/${DIS}_merge${NUM} | sed 's/NA/0/g' | sed 's/ /\t/g' > ${TEMP}/${DIS}_mutation_carriers_all.txt
  bin/transpose_perl.pl ${TEMP}/${DIS}_genes > ${TEMP}/${DIS}_genes.transpose
  cat ${TEMP}/${DIS}_mutation_carriers_all.txt ${TEMP}/${DIS}_genes.transpose > ${TEMP}/${DIS}_all_mutation_carriers.txt
  echo ""
  echo "Adding columns for MRI overlap"
  echo "ID LV" > ${DIR}/ukb_MRI_LV_ID.txt
  echo "ID RV" > ${DIR}/ukb_MRI_RV_ID.txt
  awk '{print $1, "yes"}' ${DIR}/ukb_MRI_LV_IID_clean.txt >> ${DIR}/ukb_MRI_LV_ID.txt
  awk '{print $1, "yes"}' ${DIR}/ukb_MRI_RV_IID_clean.txt >> ${DIR}/ukb_MRI_RV_ID.txt

  bin/merge_tables.pl --file1 ${DIR}/ukb_MRI_LV_ID.txt --file2 ${TEMP}/${DIS}_all_mutation_carriers.txt --index ID > ${TEMP}/${DIS}_lv
  bin/merge_tables.pl --file1 ${DIR}/ukb_MRI_RV_ID.txt --file2 ${TEMP}/${DIS}_lv --index ID > ${TEMP}/${DIS}_SNPs_MRI.txt
  echo ""
  echo "Summarizing"
  tail -n +2 ${TEMP}/${DIS}_genes | sort -u | awk '$1 != "NA" {print}' > ${DIR}/${DIS}_ExtractIID_genes.txt
  echo "" >> ${TEMP}/${DIS}_log
  echo "Genes for which there are carriers with mutations:" >> ${TEMP}/${DIS}_log
  cat ${DIR}/${DIS}_ExtractIID_genes.txt >> ${TEMP}/${DIS}_log
  echo "" >> ${TEMP}/${DIS}_log
  echo "$(awk '{print $NF}' ${TEMP}/${DIS}_SNPs_MRI.txt | grep -c yes) carriers have RV data" >> ${TEMP}/${DIS}_log
  echo "$(awk '{print $((NF-1))}' ${TEMP}/${DIS}_SNPs_MRI.txt | grep -c yes) carriers have LV data" >> ${TEMP}/${DIS}_log
  echo "" >> ${TEMP}/${DIS}_log

  Rscript --vanilla ${ROOT}/src/Extract_IID_WES_UKB.R ${TEMP}/${DIS}_SNPs_MRI.txt ${DIR}/${DIS}_ExtractIID_genes.txt ${TEMP}/${DIS}_genes.txt ${TEMP}/${DIS}_IIDs_genes_variants.txt

  cat ${TEMP}/${DIS}_log >> ${LOG}

  echo "Thanks for using this script!"
fi
