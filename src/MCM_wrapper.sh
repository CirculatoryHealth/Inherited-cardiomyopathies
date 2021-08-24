#!/bin/bash
#
#SBATCH --job-name MCM_wrapper                                            		  # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/MCM_wrapper.log              # log file of this script
#SBATCH --mem=1G                                                          		  # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                                   # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                                        # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


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
  echoerror "- Argument #1   --  Root path, from which this script will be executed, could be '/hpc/dhl_ec/mvanvugt/UKBB'"
  echoerror "- Argument #2   --  OPTIONAL: File name and path of/to the overlapping indels, could be '/hpc/dhl_ec/mvanvugt/UKBB/indels.txt' OR 'no' when this file still has to be created"
	echoerror ""
	echoerror "An example command would be: MCM_wrapper.sh [arg1: /hpc/dhl_ec/mvanvugt/UKBB] [arg3: /hpc/dhl_ec/mvanvugt/UKBB/indels.txt]."
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  	# The wrong arguments are passed, so we'll exit the script now!
  	exit 1
}

if [[ $# -lt 3 ]]; then
  echo "Error, number of arguments found "$#"."
  script_arguments_error "You must supply [3] correct arguments when running this script"

else

  ROOT="$1"
  cd ${ROOT}
  INDEL=${3:-no}
  DATE=$(date +'%d-%m-%Y')

  ### TOOLS
  SCRIPT="${ROOT}/src"
  PLINK="/hpc/local/CentOS7/dhl_ec/software/plink_v1.9"
  OVERLAP="bin/overlap.pl"
  MERGE="bin/merge_tables.pl"
  TRANSPOSE="bin/transpose_perl.pl"

  ### Directories
  UKB_200K_WES="/hpc/ukbiobank/WES_200K_2020"
  UKB_200K_WES_FAM='/hpc/dhl_ec/data/ukbiobank/WES_200K_2020'
  MRI="/hpc/dhl_ec/mvanvugt/UKBB/MRI/CMR_IDs.txt"

  ### Prepare UKB-phenotype file
  INPUT_TAB_FILE="/hpc/dhl_ec/data/ukbiobank/phenotypic/ukb44641.tab"
  INPUT_FIELD_ID="${ROOT}/data/raw/MendCM_phenotype.tsv"
  OUTPUT="${ROOT}/data/raw"
  OUTPUT_FILE_NAME="MCM"

  if [[ ! -d ${ROOT}/logs/ ]]; then
    mkdir -v ${ROOT}/logs/
  fi
  LOGS="${ROOT}/logs"
  if [[ ! -d ${ROOT}/data/temp/temp ]]; then
    mkdir -v ${ROOT}/data/temp/temp
  fi
  TEMP="${ROOT}/data/temp/temp"

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
  echo "Prefix final output: _________________________ [ ${FINAL} ]"
  echo ""
  echo "The following files will be used:"
  echo "UKB-phenotype file: __________________________ [ ${INPUT_TAB_FILE} ]"
  echo "UKB-field selection file: ____________________ [ ${INPUT_FIELD_ID} ]"
  echo "List with overlapping indels: ________________ [ ${INDEL} ]"
  echo "MRI ID list: _________________________________ [ ${MRI} ]"
  echo ""
  echo "Other options/input:"
  echo "Prefix of output: ____________________________ [ ${OUTPUT_FILE_NAME} ]"
  echo ""
  LOG="${LOGS}/${DATE}_MCM.log"
  echo "Log file created: ${LOG}"

  if [[ ! -e ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.tab ]]; then
    echo "Extracting the desired phenotypes from the full UKB-phenotype file"
    DEP1=$(sbatch --output="${LOGS}/UKBioPick.log" ${SCRIPT}/run_UKBioPick.sh ${INPUT_TAB_FILE} ${INPUT_FIELD_ID} ${OUTPUT} ${OUTPUT_FILE_NAME} --time=01:00:00 --mem 100G | sed 's/Submitted batch job //')
    echo ""
  else
    echoerror "Phenotype file ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.tab already exists. Continuing with this file. Hope that's ok!"
  fi

  echo "Compiling list of SNPs and from ClinVar and VKGL per phenotype"
  echo "Extracting participants carrying (likely) pathogenic mutations"
  DEP2a=$(sbatch --output="${LOGS}/ACM_prep_SNPs.log" ${SCRIPT}/prep_SNPs.sh ACM data/raw/ACM_Clinvar_result_LP.txt data/raw/ACM_VKGL.txt --time=30:00 | sed 's/Submitted batch job //')
  DEP3a=$(sbatch --dependency=afterok:${DEP2a} --output="${LOGS}/ACM_extract_IID.log" ${SCRIPT}/extract_IID.sh ACM ${INDEL} ${LOG} --time=01:00:00 --mem 40G | sed 's/Submitted batch job //')
  DEP2d=$(sbatch --dependency=afterok:${DEP3a} --output="${LOGS}/DCM_prep_SNPs.log" ${SCRIPT}/prep_SNPs.sh DCM data/raw/DCM_Clinvar_result_LP.txt data/raw/DCM_VKGL.txt --time=30:00 | sed 's/Submitted batch job //')
  DEP3d=$(sbatch --dependency=afterok:${DEP2d} --output="${LOGS}/DCM_extract_IID.log" ${SCRIPT}/extract_IID.sh DCM ${INDEL} ${LOG} --time=01:00:00 --mem 40G | sed 's/Submitted batch job //')
  DEP2h=$(sbatch --dependency=afterok:${DEP3d} --output="${LOGS}/HCM_prep_SNPs.log" ${SCRIPT}/prep_SNPs.sh HCM data/raw/HCM_Clinvar_result_LP.txt data/raw/HCM_VKGL.txt --time=30:00 | sed 's/Submitted batch job //')
  DEP3h=$(sbatch --dependency=afterok:${DEP2h} --output="${LOGS}/HCM_extract_IID.log" ${SCRIPT}/extract_IID.sh HCM ${INDEL} ${LOG} --time=01:00:00 --mem 40G | sed 's/Submitted batch job //')
  echo ""
  if [[ ${INDEL} == "no" ]]; then
    echobold "No file with overlapping indels provided, exiting now to make it manually or provide it"
    exit 0
  fi

  echo ""
  echo "Preparing full phenotype file"
  DEP5a=$(sbatch --dependency=afterok:${DEP3a} --output="${LOGS}/ACM_phenish.log" ${SCRIPT}/phenish.sh ACM full.tsv --time=01:30:00 --mem 100G | sed 's/Submitted batch job //')
  DEP5d=$(sbatch --dependency=afterok:${DEP3d} --output="${LOGS}/DCM_phenish.log" ${SCRIPT}/phenish.sh DCM full.tsv --time=01:30:00 --mem 100G | sed 's/Submitted batch job //')
  DEP5h=$(sbatch --dependency=afterok:${DEP3h} --output="${LOGS}/HCM_phenish.log" ${SCRIPT}/phenish.sh HCM full.tsv --time=01:30:00 --mem 100G | sed 's/Submitted batch job //')
  echo ""
  echo "Getting randomly matched IDs for the controls"
  DEP6=$(sbatch --dependency=afterok:${DEP5a},${DEP5d},${DEP5h} --output="${LOGS}/match_controls.Rlog" --wrap="Rscript --vanilla Match_controls.R data/raw _full.txt results/figures/UKB_MCM_Summary.pptx" --time=30:00)
  echo ""
  echo "Combining and cleaning up phenotypes"
  DEP7=$(sbatch --dependency=afterok:${DEP6} --output="${LOGS}/clean_pheno.log" ${SCRIPT}/clean_pheno.sh)


fi
