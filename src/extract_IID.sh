#!/bin/bash

#SBATCH --job-name extract_IID                                     # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/extract_IID.log # log file of this script
#SBATCH --mem=10G                                                  # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                      # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                           # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


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
    echoerror "- Argument #1 -- Disease to perform analysis for, could be 'DCM'"
    echoerror "- Argument #2 -- Directory to and name of configuration file, could be '/hpc/dhl_ec/mvanvugt/UKBB/config.config'"
    echoerror "- Argument #3 -- File name and path of/to the log file, could be '/hpc/dhl_ec/mvanvugt/UKBB/log.txt'"
    echoerror ""
    echoerror "An example command would be: extract_IID.sh [arg1: DCM] [arg2: /hpc/dhl_ec/mvanvugt/UKBB/config.config] [arg3: /hpc/dhl_ec/mvanvugt/UKBB/log.txt]."
    echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    # The wrong arguments are passed, so we'll exit the script now!
    exit 1
}

if [[ $# -lt 3 ]]; then
    echo "Error, number of arguments found "$#"."
    script_arguments_error "You must supply [3] correct arguments when running this script"

else

    ### Set arguments and variables
    DIS="$1"      
    CONFIG="$2"
    LOG="$3"

    # Source config and set some variables
    source ${CONFIG}
    DIR=${TEMP}
    INDEL=${RAW}/${DIS}_indels.txt

    # Make and define temporary directory
    if [[ ! -d ${DIR}/${DIS}_temp ]]; then
        mkdir -v ${DIR}/${DIS}_temp
    fi
    TEMP="${DIR}/${DIS}_temp"

    # Report on variables
    echo ""
    echo "Settings for this script:"
    echo "Disease: ____________________ [ ${DIS} ]"
    echo "Directory: __________________ [ ${DIR} ]"
    echo "Indel file: _________________ [ ${INDEL} ]"
    echo "Temporary directory: ________ [ ${TEMP} ]"
    echo "Root directory: _____________ [ ${ROOT} ]"
    echo ""

    ### STEP 2 ###
    ### overlap samples betwen 200K WES and LV/RV
    ## get IID of LV/RV MRI after QC
    tail -n +2 ${LV} | awk '{print $1,$1}' > ${TEMP}/ukb_MRI_LV_IID_clean.txt
    tail -n +2 ${RV} | awk '{print $1,$1}' > ${TEMP}/ukb_MRI_RV_IID_clean.txt

    # Adding indels to SNP-lists
    echobold "Including indels"
    echo -e "$(sort -u ${INDEL} | wc -l | cut -d" " -f1) indels included" >> ${TEMP}/${DIS}_log
    awk '{print $1}' ${INDEL} >> ${DIR}/${DIS}_overlap_LP_WES_SNPs.txt

    ## Allele freq. using Plink
    echobold "Calculating allele frequencies"
    rm ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.allele_frq
    for CHR in $(seq 1 22); do
        echo ${CHR}
        ${PLINK} --bed ${UKB_200K_WES}/${FIL_NAM}${CHR}${ADD}.bed \
            --bim  ${UKB_200K_WES}/${FIL_NAM}${CHR}${ADD}.bim \
            --fam  ${UKB_200K_WES_FAM}/${FIL_NAM}${CHR}${ADD}.fam \
            --extract ${DIR}/${DIS}_overlap_LP_WES_SNPs.txt \
            --freq --out ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR} --silent
        tail -n +2 ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.frq >> ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.allele_frq
        echo ""
    done

    # Remove SNPs without carriers and report on number of SNPs with carriers
    echobold "Including SNPs with carriers only"
    awk '$5 != 0 {print $0}' ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.allele_frq > ${DIR}/${DIS}_WES_MRI_UKB_chrALL.allele_frq_filtered
    awk '{print $2}' ${DIR}/${DIS}_WES_MRI_UKB_chrALL.allele_frq_filtered > ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt
    echo -e "$(sort -u ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt | wc -l | cut -d" " -f1) SNPs with carriers in the UK Biobank" >> ${TEMP}/${DIS}_log

    # Remove SNPs that have overlap between HCM/DCM/ACM and are not annotated well
    # Extract all to be included SNPs from VKGL and check in Clinvar. Annotated for ${DIS} of cardiomyopathy? Include, else exclude -- add column 2 to data/raw/${DIS}_remove_overlap.txt
    echobold "Removing SNPs specified in ${RAW}/${DIS}_remove_overlap.txt"
    ${OVERLAP} ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt 1 ${TEMP}/${DIS}_VKGL_positionID 1 > ${TEMP}/${DIS}_VKGL_check
    if [[ -e ${RAW}/${DIS}_remove_overlap.txt ]]; then
        echo -e "$(sort -u ${RAW}/${DIS}_remove_overlap.txt | wc -l | cut -d" " -f1) SNPs removed due to ambiguous annotation" >> ${TEMP}/${DIS}_log
        ${OVERLAP} ${RAW}/${DIS}_remove_overlap.txt 1 ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt 1 -v > ${TEMP}/${DIS}_temp
        mv ${TEMP}/${DIS}_temp ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt
    else
        echoerror "${RAW}/${DIS}_remove_overlap.txt does not exist, not removing overlapping SNPs"
    fi

    # Report on final number of SNPs
    echo -e "$(sort -u ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt | wc -l | cut -d" " -f1) SNPs to extract UKB-individuals for" >> ${TEMP}/${DIS}_log
    echo "" >> ${TEMP}/${DIS}_log
    # Create VCF-files for final SNP-list
    echobold "Making vcf-files with only desired snps"
    for CHR in $(seq 1 22); do
        rm ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.vcf
        echo ${CHR}
        ${PLINK} --bed ${UKB_200K_WES}/${FIL_NAM}${CHR}${ADD}.bed \
            --bim  ${UKB_200K_WES}/${FIL_NAM}${CHR}${ADD}.bim \
            --fam  ${UKB_200K_WES_FAM}/${FIL_NAM}${CHR}${ADD}.fam \
            --extract ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt \
            --recode vcf-iid --silent \
            --keep-allele-order --out ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}
        echo ""
    done

    echo "Making one vcf-file"
    # Save the number of a chromosome definitely in this dataset
    HEADER_CHR_2=$(head -1 ${DIR}/${DIS}_overlap_LP_WES_SNPs_updated.txt | cut -d':' -f1)
    # Create header using one of the vcf-files
    head -7 ${TEMP}/${DIS}_WES_MRI_UKB_chr${HEADER_CHR_2}.vcf > ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.vcf
    # Combine all vcf-files
    for CHR in $(seq 1 22); do
        tail -n +8 ${TEMP}/${DIS}_WES_MRI_UKB_chr${CHR}.vcf >> ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.vcf
    done

    echo "Making list of genes per SNP"
    # Remove meta-information
    awk '$1 !~ /^##/' ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.vcf > ${TEMP}/${DIS}_chrALL.vcf
    # Create list of SNPs with corresponding genes
    echo "ID Gene" > ${TEMP}/${DIS}_sng_temp
    tail -n +2 ${TEMP}/${DIS}_LP_positionID | awk '{print $2, $4}' >> ${TEMP}/${DIS}_sng_temp
    cat ${INDEL} >> ${TEMP}/${DIS}_sng_temp
    # Add Gene information to vcf-file
    ${OVERLAP} ${TEMP}/${DIS}_chrALL.vcf 3 ${TEMP}/${DIS}_sng_temp 1 | sort -ur > ${TEMP}/snp_gen
    ${MERGE} --file1 ${TEMP}/${DIS}_chrALL.vcf --file2 ${TEMP}/snp_gen --index ID | sed 's/ /\t/g' | cut -f1,2,11- > ${TEMP}/${DIS}_gene.vcf
    # Isolate Genes
    cut -f2 ${TEMP}/${DIS}_gene.vcf > ${TEMP}/${DIS}_genes
    echo ""
    echo "Getting SNP and gene per individual"
    # Get SNP and gene information
    tail -n +2 ${TEMP}/${DIS}_gene.vcf | cut -f1,2 > ${TEMP}/${DIS}_SNPs.txt
    # Initiate file if necessary
    if [[ ! -e ${ROOT}/data/processed/All_SNP_IID.txt ]]; then
        echo "f.eid SNP Gene" > ${ROOT}/data/processed/All_SNP_IID.txt
    fi

    # Read ${TEMP}/${DIS}_SNPs.txt line by line
    while IFS= read -r line; do

        # Save content of line
        entries=($line)
        SNP=$(echo "${entries[0]}")
        GENE=$(echo "${entries[1]}")

        # Per SNP get carriers?
        echo "${SNP}"
        NUM=`grep ${SNP} ${TEMP}/${DIS}_gene.vcf | ${TRANSPOSE} | grep -n 1 | tail -n +2 | sed 's/:/\t/g' | cut -f1 | ${TRANSPOSE} | sed 's/\t/,/g'`
        head -1 ${TEMP}/${DIS}_gene.vcf | cut -f${NUM} | ${TRANSPOSE} | awk -v snp="${SNP}" -v gen="${GENE}" '{print $1, snp, gen}' >> ${ROOT}/data/processed/All_SNP_IID.txt

    done < "${TEMP}/${DIS}_SNPs.txt"
    echo ""
    echo "Transposing vcf-file so header are SNPs and rows are individuals"
    cat ${TEMP}/${DIS}_gene.vcf | ${TRANSPOSE} > ${TEMP}/${DIS}_WES_MRI_UKB_chrALL.transpose

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
            ${MERGE} --file1 ${TEMP}/${DIS}_SNP${SNP}_ID.txt --file2 ${TEMP}/${DIS}_LP_allIDs.txt --index ID > ${TEMP}/${DIS}_merge${SNP}
        else
            MIN=$((${SNP}-1))
            ${MERGE} --file1 ${TEMP}/${DIS}_SNP${SNP}_ID.txt --file2 ${TEMP}/${DIS}_merge${MIN} --index ID > ${TEMP}/${DIS}_merge${SNP}
        fi
    done

    # Clean up
    sed 's/0\///g' ${TEMP}/${DIS}_merge${NUM} | sed 's/ NA/ 0/g' | sed 's/ /\t/g' > ${TEMP}/${DIS}_mutation_carriers_all.txt
    echo ""
    echo "Adding columns for CMR availability"
    echo "ID LV" > ${DIR}/ukb_MRI_LV_ID.txt
    echo "ID RV" > ${DIR}/ukb_MRI_RV_ID.txt
    awk '{print $1, "yes"}' ${TEMP}/ukb_MRI_LV_IID_clean.txt >> ${DIR}/ukb_MRI_LV_ID.txt
    awk '{print $1, "yes"}' ${TEMP}/ukb_MRI_RV_IID_clean.txt >> ${DIR}/ukb_MRI_RV_ID.txt

    ${MERGE} --file1 ${DIR}/ukb_MRI_LV_ID.txt --file2 ${TEMP}/${DIS}_mutation_carriers_all.txt --index ID > ${TEMP}/${DIS}_lv
    ${MERGE} --file1 ${DIR}/ukb_MRI_RV_ID.txt --file2 ${TEMP}/${DIS}_lv --index ID > ${TEMP}/${DIS}_SNPs_MRI.txt
    echo ""
    echo "Summarizing"
    tail -n +2 ${TEMP}/${DIS}_genes | sort -u | awk '$1 != "NA" {print}' > ${DIR}/${DIS}_ExtractIID_genes.txt
    echo "" >> ${TEMP}/${DIS}_log
    echo "Genes for which there are carriers with mutations:" >> ${TEMP}/${DIS}_log
    cat ${DIR}/${DIS}_ExtractIID_genes.txt >> ${TEMP}/${DIS}_log
    echo "" >> ${TEMP}/${DIS}_log
    echo -e "$(tail -n +3 ${TEMP}/${DIS}_SNPs_MRI.txt | sort -u | wc -l | cut -d" " -f1) individuals carrying these ${DIS}-associated mutations included" >> ${TEMP}/${DIS}_log
    echo "$(awk '{print $NF}' ${TEMP}/${DIS}_SNPs_MRI.txt | grep -c yes) carriers have RV data" >> ${TEMP}/${DIS}_log
    echo "$(awk '{print $((NF-1))}' ${TEMP}/${DIS}_SNPs_MRI.txt | grep -c yes) carriers have LV data" >> ${TEMP}/${DIS}_log
    echo "" >> ${TEMP}/${DIS}_log

    # Summarize variant-level data on gene-level
    ${ROOT}/src/Extract_IID_WES_UKB.R ${TEMP}/${DIS}_SNPs_MRI.txt ${DIR}/${DIS}_ExtractIID_genes.txt ${TEMP}/${DIS}_genes.txt ${TEMP}/${DIS}_IIDs_genes_variants.txt

    cat ${TEMP}/${DIS}_log >> ${LOG}
    rm ${TEMP}/${DIS}_log
    echo "Thanks for using this script!"
fi
