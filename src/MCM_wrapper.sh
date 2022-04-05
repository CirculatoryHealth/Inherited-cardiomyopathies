#!/bin/bash
#
#SBATCH --job-name MCM_wrapper                                     # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/MCM_wrapper.log # log file of this script
#SBATCH --mem=1G                                                   # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                      # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                           # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                           Wrapper to perform analyses on HCM/DCM/ACM pathogenic mutation carriers                            "
echo "                                                    version 1.0 (20210714)                                                    "
echo ""
echo "* Written by      : Marion van Vugt"
echo "* E-mail          : m.vanvugt-2@umcutrecht.nl"
echo "* Last update     : 2021-07-14"
echo "* Version         : MCM_wrapper_1.0"
echo ""
echo "* Description     : This script will execute scripts that perform analyses on HCM/DCM/ACM pathogenic mutation carriers     "
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
    echoerror "- Argument #1 -- Configuration file, could be '/hpc/dhl_ec/mvanvugt/UKBB/config/config'"
    echoerror ""
    echoerror "An example command would be: MCM_wrapper.sh [arg1: /hpc/dhl_ec/mvanvugt/UKBB/config/config]."
    echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    # The wrong arguments are passed, so we'll exit the script now!
    exit 1
}

if [[ $# -lt 1 ]]; then
    echo "Error, number of arguments found "$#"."
    script_arguments_error "You must supply [1] correct arguments when running this script"

else

    CONFIG="$1"
    source ${CONFIG}
    cd ${ROOT}
    DATE=$(date +'%d-%m-%Y')

    ### Make and define some directories
    OUTPUT="${RAW}"
    if [[ ! -d ${ROOT}/logs/ ]]; then
        mkdir -v ${ROOT}/logs/
    fi
    LOGS="${ROOT}/logs"
    if [[ ! -d ${ROOT}/data/temp/temp ]]; then
        mkdir -v ${ROOT}/data/temp/temp
    fi
    TEMP="${ROOT}/data/temp/temp"

    ### Report on all variables
    echo ""
    echo "Setting directories:"
    echo "Script directory: ____________________________ [ ${SCRIPT} ]"
    echo "Root directory: ______________________________ [ ${ROOT} ]"
    echo "UKB WES genetic data: ________________________ [ ${UKB_200K_WES} ]"
    echo "UKB WES fam-directory: _______________________ [ ${UKB_200K_WES_FAM} ]"
    echo "Output directory: ____________________________ [ ${OUTPUT} ]"
    echo "Log-file directory: __________________________ [ ${LOGS} ]"
    echo "Temporary directory: _________________________ [ ${TEMP} ]"
    echo ""
    echo "The following software will be used:"
    echo "Plink: _______________________________________ [ ${PLINK} ]"
    echo "Overlap tool: ________________________________ [ ${OVERLAP} ]"
    echo "Merge tool: __________________________________ [ ${MERGE} ]"
    echo "Transpose tool: ______________________________ [ ${TRANSPOSE} ]"
    echo ""
    echo "The following files will be used:"
    echo "UKB-phenotype file: __________________________ [ ${INPUT_TAB_FILE} ]"
    echo "UKB-field selection file: ____________________ [ ${INPUT_FIELD_ID} ]"
    echo "CMR full data: _______________________________ [ ${FULL_CMR} ]"
    echo ""
    echo "Other options/input:"
    echo "Prefix of output: ____________________________ [ ${OUTPUT_FILE_NAME} ]"
    echo ""
    LOG="${LOGS}/${DATE}_MCM.log"
    echo "Log file created: ${LOG}"

    # If necessary, extract phenotypes from UKB
    if [[ ! -e ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.tab ]]; then
        echo "Extracting the desired phenotypes from the full UKB-phenotype file"
        DEP1=$(sbatch --output="${LOGS}/UKBioPick.log" ${ROOT}/bin/UKBioPick.sh ${INPUT_TAB_FILE} ${INPUT_FIELD_ID} ${OUTPUT} ${OUTPUT_FILE_NAME} --time=01:00:00 --mem 100G | sed 's/Submitted batch job //')
        echo ""
    else
        echoerror "Phenotype file ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.tab already exists. Continuing with this file. Hope that's ok!"
    fi

    # Create CMR ID list
    awk '{print $1}' ${FULL_CMR} > ${ROOT}/data/temp/CMR_IDs.txt

    # Initiate jobs per CM to compile list of variants and participants using the downloaded input from ClinVar and VKGL
    echo "Compiling list of SNPs and from ClinVar and VKGL per phenotype"
    echo "Extracting participants carrying (likely) pathogenic mutations"
    DEP2a=$(sbatch --time=30:00 --output="${LOGS}/ACM_prep_SNPs.log" ${SCRIPT}/prep_SNPs.sh ACM ${CONFIG} ${LOG} | sed 's/Submitted batch job //')
    DEP3a=$(sbatch --time=01:00:00 --mem 40G --dependency=afterok:${DEP2a} --output="${LOGS}/ACM_extract_IID.log" ${SCRIPT}/extract_IID.sh ACM ${CONFIG} ${LOG} | sed 's/Submitted batch job //')
    DEP2d=$(sbatch --time=30:00 --dependency=afterok:${DEP3a} --output="${LOGS}/DCM_prep_SNPs.log" ${SCRIPT}/prep_SNPs.sh DCM ${CONFIG} ${LOG} | sed 's/Submitted batch job //')
    DEP3d=$(sbatch --time=01:00:00 --mem 40G --dependency=afterok:${DEP2d} --output="${LOGS}/DCM_extract_IID.log" ${SCRIPT}/extract_IID.sh DCM ${CONFIG} ${LOG}| sed 's/Submitted batch job //')
    DEP2h=$(sbatch --time=30:00 --dependency=afterok:${DEP3d} --output="${LOGS}/HCM_prep_SNPs.log" ${SCRIPT}/prep_SNPs.sh HCM ${CONFIG} ${LOG} | sed 's/Submitted batch job //')
    DEP3h=$(sbatch --time=01:00:00 --mem 40G --dependency=afterok:${DEP2h} --output="${LOGS}/HCM_extract_IID.log" ${SCRIPT}/extract_IID.sh HCM ${CONFIG} ${LOG} | sed 's/Submitted batch job //')
    echo ""

    # Initiate job per CM to create the phenotype file
    echo ""
    echo "Preparing full phenotype file"
    DEP5a=$(sbatch --time=01:30:00 --mem 100G --dependency=afterok:${DEP3a} --output="${LOGS}/ACM_phenish.log" ${SCRIPT}/phenish.sh ACM ${CONFIG} full.tsv | sed 's/Submitted batch job //')
    DEP5d=$(sbatch --time=01:30:00 --mem 100G --dependency=afterok:${DEP3d} --output="${LOGS}/DCM_phenish.log" ${SCRIPT}/phenish.sh DCM ${CONFIG} full.tsv | sed 's/Submitted batch job //')
    DEP5h=$(sbatch --time=01:30:00 --mem 100G --dependency=afterok:${DEP3h} --output="${LOGS}/HCM_phenish.log" ${SCRIPT}/phenish.sh HCM ${CONFIG} full.tsv | sed 's/Submitted batch job //')

    # Extract phenotypes required to match controls to variant-carriers
    ${ROOT}/bin/UKBioPick.sh ${INPUT_TAB_FILE} ${CONTROL_FIELDS} ${ROOT}/data/temp MCM_WES
    # solve delimiter problem
    sed 's/ /_/g' ${ROOT}/data/temp/MCM_WES_ukb_phenotypes.tab > ${ROOT}/data/temp/MCM_WES_pheno.tsv
    # Include controls
    echo ""
    echo "Getting randomly matched IDs for the controls"
    DEP6=$(sbatch --time=30:00 --dependency=afterok:${DEP5a},${DEP5d},${DEP5h} --output="${LOGS}/match_controls.Rlog" --wrap="${SCRIPT}/Match_controls.R ${OUTPUT} _full.tsv results/figures ${MUL}")
    echo ""

    # Create full and clean phenotype file
    echo "Combining and cleaning up phenotypes"
    mv ${ROOT}/data/temp/ACM_temp/WES_MRI_MCM_phenotypes.txt ${ROOT}/data/temp
    DEP7=$(sbatch --dependency=afterok:${DEP6} --output="${LOGS}/clean_pheno.log" ${SCRIPT}/clean_pheno.sh ${CONFIG})


fi
