#!/bin/bash
#
#SBATCH --job-name clean_pheno                                     # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/clean_pheno.log # log file of this script
#SBATCH --mem=10G                                                  # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                      # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                           # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                                       Script to combine and clean up the phenotype file                                      "
echo "                                                    version 1.0 (20210714)                                                    "
echo ""
echo "* Written by      : Marion van Vugt"
echo "* E-mail          : m.vanvugt-2@umcutrecht.nl"
echo "* Last update     : 2021-07-14"
echo "* Version         : clean_pheno_1.0"
echo ""
echo "* Description     : This script will combine and clean up the phenotype file      "
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
    echoerror "An example command would be: clean_pheno.sh [arg1: /hpc/dhl_ec/mvanvugt/UKBB/config/config]."
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

    echo "Making phenotype file for the controls"
    ${MERGE} --file1 ${TEMP}/WES_MRI_MCM_phenotypes.txt --file2 ${RAW}/Control_IDs.tsv --index f.eid | sed 's/ /\t/g' > ${RAW}/Controls_full.tsv
    echo ""
    echo "Combining the phenotype files"
    ${ROOT}/src/Combine_pheno.R ${RAW} _full.tsv MCM_raw_full
    echo ""
    echo "Summarizing genetic information"
    ${ROOT}/src/MCM_gene_summary.R ${RAW}/MCM_raw_full.rds ${TEMP} MCM
    ${MERGE} --file1 ${ROOT}/data/processed/All_SNP_IID.txt --file2 ${TEMP}/MCM_gene_summary.tsv --index f.eid | sort -ur | sed 's/ /\t/g' > ${ROOT}/data/processed/MCM_gene_summary.tsv
    echo ""
    echo "Cleaning up the phenotype file"
    ${ROOT}/src/MCM_pheno_clean.R ${RAW}/MCM_raw_full.rds ${RAW} ${ROOT}/data/processed/ MCM
    echo ""
    echo "Merging the cleaned files"
    ${MERGE} --file1 ${ROOT}/data/processed/MCM_gene_summary.tsv --file2 ${ROOT}/data/processed/MCM_cleaned.tsv --index f.eid | sed 's/ /\t/g' > ${ROOT}/data/processed/MCM_final_pheno.tsv

fi
