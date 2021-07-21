#!/bin/bash
#
#SBATCH --job-name prep_SNPs                                              		  # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/prep_SNPs.log                # log file of this script
#SBATCH --mem=10G                                                          		  # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                                   # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                                        # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                                     Script to prepare the list of SNPs associated with CM                                    "
echo "                                                    version 1.0 (20210714)                                                    "
echo ""
echo "* Written by      : Marion van Vugt"
echo "* E-mail          : m.vanvugt-2@umcutrecht.nl"
echo "* Last update     : 2021-07-14"
echo "* Version         : prep_SNPs_1.0"
echo ""
echo "* Description     : This script will prepare the list of SNPs associated with CM     "
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
  echoerror "- Argument #2   --  Name of and path to file with the ClinVar results, could be 'data/raw/DCM_Clinvar_result_LP.txt'"
  echoerror "- Argument #3   --  Name of and path to file with the VKGL results, could be 'DCM_VKGL.txt'"
  echoerror "- Argument #4   --  Name of and path to the log-file, could be 'data/raw/log.txt'"
	echoerror ""
	echoerror "An example command would be: prep_SNPs.sh [arg1: DCM] [arg2: DCM_Clinvar_result_LP.txt] [arg3: DCM_VKGL.txt] [arg4: data/raw/log.txt]."
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  	# The wrong arguments are passed, so we'll exit the script now!
  	exit 1
}

if [[ $# -lt 4 ]]; then
  echo "Error, number of arguments found "$#"."
  script_arguments_error "You must supply [4] correct arguments when running this script"

else

  DIS="$1"      ### Disease (DCM/HCM/ACM)
  CV="$2"
  VKGL="$3"
  LOG="$4"
  DIR="data/temp"
  OVERLAP="bin/overlap.pl"

  if [[ ! -d ${DIR}/${DIS}_temp ]]; then
    mkdir -v ${DIR}/${DIS}_temp
  else
    echo "Temporary folder already exists, emptying it now"
    rm ${DIR}/${DIS}_temp/*
  fi
  TEMP="${DIR}/${DIS}_temp"

  echo ""
  echo "Settings for this script:"
  echo "Disease: ____________________ [ ${DIS} ]"
  echo "Directory: __________________ [ ${DIR} ]"
  echo "Clinvar file: _______________ [ ${CV} ]"
  echo "VKGL file: __________________ [ ${VKGL} ]"
  echo "Temporary directory: ________ [ ${TEMP} ]"
  echo "Log-file: ___________________ [ ${LOG} ]"
  echo ""

  ### STEP 1 ###

  ## In ClinVar search ""Arrhythmogenic right ventricular cardiomyopathy"[dis] OR "Arrhythmogenic right ventricular dysplasia[dis] OR "Arrhythmogenic cardiomyopathy"[dis]" ## put all alternative name of the disease
  ## Click on Pathogenic and Likely Pathogenic
  ## Download file with these features; Format:Tabular (text), Sort by:Location
  ## Upload file to HPC and rename it with ${DIS}_clinvar_result_LP.txt

  echo "=======================================================================" > ${TEMP}/${DIS}_log
  echo "                         Extraction ${DIS} SNPs                        " >> ${TEMP}/${DIS}_log
  echo "=======================================================================" >> ${TEMP}/${DIS}_log
  echo "" >> ${TEMP}/${DIS}_log

  echo "Extract SNPs from the databases"
  echo -e "$(wc -l ${CV} | cut -d" " -f1) SNPs extracted from the ClinVar database" >> ${TEMP}/${DIS}_log
  echo -e "$(wc -l ${VKGL} | cut -d" " -f1) SNPs extracted from the VKGL database" >> ${TEMP}/${DIS}_log

  echo "ID Position_38 Position_37 Gene" > ${TEMP}/${DIS}_LP_positionID
  tail -n +2 ${CV} | sed 's/ /_/g' | awk -F"\t" '{print $2, $8, $9, $10, $11, $15}' | sed 's/:/ /g' | awk '{print $4 ":" $7 ":" $8 ":" $9, $4 ":" $5 ":" $8 ":" $9, $2 ":" $3 ":" $8 ":" $9, $1}' | sed 's/_-_\([0-9]\+\):/:/g' > ${TEMP}/${DIS}_LP_positionID_temp
  tail -n +2 ${VKGL} | awk '{print $3 ":" $9 ":" $10, $3 ":" $9 ":" $10, $2 ":" $9 ":" $10, $4}' >> ${TEMP}/${DIS}_LP_positionID_temp
  echo ""
  # Filter for certain genes
  awk -v cm=${DIS} '$2 == cm {print $1}' ${DIR}/CM_genes.txt > ${TEMP}/${DIS}_genes.txt
  echo "Checking for gene: "
  while IFS= read -r line; do

    IFS=" " read -ra entries <<< "$line"
    GENE=$(echo "${entries[0]}")

    echo "${GENE}"
    # Only add the variants that are in certain genes
    grep -w ${GENE} ${TEMP}/${DIS}_LP_positionID_temp >> ${TEMP}/${DIS}_LP_positionID

  done < "${TEMP}/${DIS}_genes.txt"
  echo -e "$(wc -l ${TEMP}/${DIS}_LP_positionID | cut -d" " -f1) SNPs remaining after filtering for ${DIS}-associated genes" >> ${TEMP}/${DIS}_log
  echo ""

  echo "Extracting SNPs present in the WES-data"
  ## Overlap with WES
  ${OVERLAP} ${TEMP}/${DIS}_LP_positionID 2 data/raw/ukb_SNP_ID_WES_chrALL.txt 2 | cut -d' ' -f1 | sort -u >  ${DIR}/${DIS}_overlap_LP_WES_SNPs.txt
  sed 's/:/ /g5' data/raw/ukb_SNP_ID_WES_chrALL.txt > data/raw/ukb_SNP_ID_WES_chrALL.txt_position # replace ":" between position and ref to space
  sed 's/:/ /g2' ${TEMP}/${DIS}_LP_positionID | awk '{print $1, $4 ":" $5}' > ${TEMP}/${DIS}_LP_positionID_position
  ${OVERLAP} ${TEMP}/${DIS}_LP_positionID_position 2 data/raw/ukb_SNP_ID_WES_chrALL.txt_position 2 | sort -u > ${TEMP}/${DIS}_overlap_LP_WES_SNPs.txt_position # get overlap by position
  echo "" >> ${TEMP}/${DIS}_log
  cat ${TEMP}/${DIS}_log >> ${LOG}
  rm ${TEMP}/${DIS}_log
  echo ""
  echo "Finished! Ciao!"

fi
