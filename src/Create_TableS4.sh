#!/bin/bash
#
#SBATCH --job-name Create_TableS4                                     # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/Create_TableS4.log # log file of this script
#SBATCH --time=1:30:00                                                # --time=[max time, e.g. 02:02:01] - this is the time you think the script will take
#SBATCH --mem=30G                                                     # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                         # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                              # you can choose: ALL=select all options; END=end of job; FAIL=abort of job

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                            Script to start creating Table S4                           "
echo "                                 version 1.0 (20220415)                                 "
echo ""
echo "* Written by      : Marion van Vugt"
echo "* E-mail          : m.vanvugt-2@umcutrecht.nl"
echo "* Last update     : 2022-04-15"
echo "* Version         : Create_TableS4_1.0"
echo ""
echo "* Description     : This script will start creating Table S4"
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
date "+DATE: %a %d/%m/%Y%nTIME: %T"
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
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
    echoerror "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echoerror "$1" # ERROR MESSAGE
    echoerror "- Argument #1 -- Path to and name of the configuration file, could be '/hpc/dhl_ec/mvanvugt/UKBB/config/config'"
    echoerror ""
    echoerror "An example command would be: Create_TableS4.sh [arg1: /hpc/dhl_ec/mvanvugt/UKBB/config/config]"
    echoerror "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    # The wrong arguments are passed, so we'll exit the script now!
    exit 1
}

if [[ $# -lt 1 ]]; then
    echo "Error, number of arguments found "$#"."
    script_arguments_error "You must supply [1] correct arguments when running this script"

else

    # Specifying and sourcing configuration file
    CONFIG="$1"
    source ${CONFIG}

    # Reporting on all set variables
    echo "Setting directories, files and options:"
    echo "Name without path first argument: _____________________ [ ${INPUT##/*/} ]"
    echo "Path without name first argument: _____________________ [ ${INPUT%/*} ]"
    echo "Root directory: _______________________________________ [ ${ROOT} ]"
    echo "Directory to temporary files: _________________________ [ ${TEMP} ]"
    echo "Directory to raw data files: __________________________ [ ${RAW} ]"
    echo "Script directory: _____________________________________ [ ${SCRIPT} ]"
    echo ""

    echobold "Make table of all SNPs with their respective numbers per CM"
    # Create full list of SNPs included
    echo "SNP" > ${TEMP}/CM_incl_SNPs.txt
    cat ${TEMP}/ACM_overlap_LP_WES_SNPs_updated.txt ${TEMP}/DCM_overlap_LP_WES_SNPs_updated.txt ${TEMP}/HCM_overlap_LP_WES_SNPs_updated.txt | sort -u >> ${TEMP}/CM_incl_SNPs.txt

    # Per CM calculate the number of carriers per SNP
    for DIS in ACM DCM HCM; do

        # Print header
        echo "SNP N_${DIS}"> ${TEMP}/overlap_SNPs_${DIS}.txt

        # Read SNP-list of this disease and count number of carriers included
        while IFS= read -r line; do
            # Search for SNP and count number of carriers
            echo "Searching for the following SNP now: ${line}"
            COL=`head -1 ${TEMP}/${DIS}_temp/${DIS}_SNPs_MRI.txt | ${TRANSPOSE} | grep -n ${line} | sed 's/:/\t/' | cut -f1`
            # Save number of carriers and SNP position to file
            echo -e "${line} $(awk -v col="${COL}" '{print $col}' ${TEMP}/${DIS}_temp/${DIS}_SNPs_MRI.txt | tail -n +2 | grep -c 1)" >> ${TEMP}/overlap_SNPs_${DIS}.txt

        done < "${TEMP}/${DIS}_overlap_LP_WES_SNPs_updated.txt"

        echobold "Add other important information"
        # Solve some delimiter problems
        sed 's/:/\t/2' ${TEMP}/CM_incl_SNPs.txt > ${TEMP}/CM_incl
        sed -i 's/ /_/g' ${RAW}/${DIS}_clinvar_result_LP.txt
        sed -i 's/ /_/g' ${RAW}/${DIS}_VKGL.txt

        # Create correct position ID for ClinVar SNPs
        cut -f 10,11,15 ${RAW}/${DIS}_clinvar_result_LP.txt | sed 's/:/\t/2' | awk '{print $1 ":" $2 ":" $4, $3 ":" $4}' > ${TEMP}/${DIS}_SNP
        # Create new file for position ID and canonical SPDI
        echo "SNP Canonical_SPDI" > ${TEMP}/${DIS}_SNP_info
        # Populate file with position ID and canonical SPDI
        ${OVERLAP} ${TEMP}/overlap_SNPs_${DIS}.txt 1 ${TEMP}/${DIS}_SNP 1 >> ${TEMP}/${DIS}_SNP_info

        # Include even more information from the ClinVar SNPs
        echo "SNP GRCh37 Gene Canonical_SPDI Origin Pathogenicity" | sed 's/ /\t/g' > ${TEMP}/${DIS}_all_SNP
        ${MERGE} --file1 ${RAW}/${DIS}_clinvar_result_LP.txt --file2 ${TEMP}/${DIS}_SNP_info --index Canonical_SPDI | awk 'NR > 1 {print $1, $10 ":" $11, $4, $2, "ClinVar", $7}' | sed 's/ /\t/g' >> ${TEMP}/${DIS}_all_SNP

        # Create correct position ID for VKGL SNP
        cut -f 2-4,9,10,14,16 ${RAW}/${DIS}_VKGL.txt | awk '{print $2 ":" $4 ":" $5, $1, $3, $6, "VKGL", $7}' > ${TEMP}/${DIS}_VKGL
        # Create list of only included SNPs
        ${OVERLAP} ${TEMP}/overlap_SNPs_${DIS}.txt 1 ${TEMP}/${DIS}_VKGL 1 > ${TEMP}/${DIS}_VKGL_info
        ${OVERLAP} ${TEMP}/${DIS}_all_SNP 1 ${TEMP}/${DIS}_VKGL_info 1 -v | sed 's/ /\t/g' | tail -n +2 >> ${TEMP}/${DIS}_all_SNP

        echoerror "At this point, please manually add the indels to this list"
        # Manually add the indels to this list using the following commands:
        # grep loc data/raw/${DIS}_clinvar_result_LP.txt | awk '{print SNP_38, $8 ":" $9, $2, $15, "ClinVar", $5}' | sed 's/ /\t/g' >> data/temp/${DIS}_all_SNP
        # grep loc data/raw/${DIS}_VKGL.txt | awk '{print SNP_38, $2, $4, $14, "VKGL", $16}' | sed 's/ /\t/g' >> data/temp/${DIS}_all_SNP
        done

    echobold "Combining the number of carriers and extra information"
    # Add the Number of individuals per CM to full SNP list
    ${MERGE} --file1 ${TEMP}/overlap_SNPs_ACM.txt --file2 ${TEMP}/CM_incl_SNPs.txt --index SNP > ${TEMP}/tempACM
    ${MERGE} --file1 ${TEMP}/overlap_SNPs_DCM.txt --file2 ${TEMP}/tempACM --index SNP > ${TEMP}/tempDCM
    ${MERGE} --file1 ${TEMP}/overlap_SNPs_HCM.txt --file2 ${TEMP}/tempDCM --index SNP | sed 's/ /\t/g' > ${TEMP}/tempHCM


    # Add the extra info to full SNP list
    cat ${TEMP}/ACM_all_SNP ${TEMP}/DCM_all_SNP ${TEMP}/HCM_all_SNP | sort -ur > ${TEMP}/All_SNP
    ${MERGE} --file1 ${TEMP}/All_SNP --file2 ${TEMP}/tempHCM --index SNP > ${TEMP}/SNP_table

    echo "Add MAF to table"
    echo "SNP MAF" > ${TEMP}/SNPs_MAF.frq
    # Combine all MAF information
    cat ${TEMP}/ACM_WES_MRI_UKB_chrALL.allele_frq_filtered ${TEMP}/DCM_WES_MRI_UKB_chrALL.allele_frq_filtered ${TEMP}/HCM_WES_MRI_UKB_chrALL.allele_frq_filtered | awk '{print $2, $5}' | sort -u >> ${TEMP}/SNPs_MAF.frq
    # Merge MAF with table
    echobold "Saving the pre-Table S4 in ${ROOT}/results/output/TableS4.tsv"
    ${MERGE} --file1 ${TEMP}/SNPs_MAF.frq --file2 ${TEMP}/SNP_table --index SNP | sed 's/ /\t/g' | sed 's/ACM/ARVC/' > ${ROOT}/results/output/TableS4.tsv

    echoerror "The rest of Table S4 was created manually"
    # Cleaning up
    rm ${TEMP}/CM_incl* ${TEMP}/overlap_SNPs_ACM.txt ${TEMP}/overlap_SNPs_DCM.txt ${TEMP}/overlap_SNPs_HCM.txt ${TEMP}/temp* ${TEMP}/*_VKGL* ${TEMP}/*_SNP_* ${TEMP}/*_SNP ${TEMP}/SNP_table ${TEMP}/SNPs_MAF.frq

    # Announce the end of the script
    echobold "Finished, thanks for using this script!"

fi

script_copyright_message


