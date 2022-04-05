#!/bin/bash
#
#SBATCH --job-name phenish                                     # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/phenish.log # log file of this script
#SBATCH --mem=10G                                              # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                  # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                       # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                                    Script to prepare the full phenotype file per phenotype                                   "
echo "                                                    version 1.0 (20210714)                                                    "
echo ""
echo "* Written by      : Marion van Vugt"
echo "* E-mail          : m.vanvugt-2@umcutrecht.nl"
echo "* Last update     : 2021-07-14"
echo "* Version         : phenish_1.0"
echo ""
echo "* Description     : This script will prepare the full phenotype file per phenotype     "
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
    echoerror "- Argument #2 -- Configuration file, could be '/hpc/dhl_ec/mvanvugt/UKBB/Inherited-cardiomyopathies/config/config'"
    echoerror "- Argument #3 -- Suffix of the output file, could be 'full.tsv'"
    echoerror ""
    echoerror "An example command would be: phenish.sh [arg1: DCM] [arg2: /hpc/dhl_ec/mvanvugt/UKBB/Inherited-cardiomyopathies/data/raw/MCM_ukb_phenotypes.tab] [arg3: full.txt]."
    echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    # The wrong arguments are passed, so we'll exit the script now!
    exit 1
}

if [[ $# -lt 3 ]]; then
    echo "Error, number of arguments found "$#"."
    script_arguments_error "You must supply [3] correct arguments when running this script"

else

    # Set arguments and source configuration file
    DIS="$1"      
    CONFIG="$2"
    SUF="$3"
    source ${CONFIG}

    # Set variables
    DIR="${TEMP}"
    OUT="${RAW}"
    TEMP="${DIR}/${DIS}_temp"
    PHENO="${RAW}/${OUTPUT_FILE_NAME}_ukb_phenotypes.tab"

    # Clean up some details of the file (delimiter and ID-name)
    sed 's/ /\t/g' ${TEMP}/${DIS}_IIDs_genes_variants.txt | sed 's/ID/f.eid/' > ${TEMP}/${DIS}_IIDs_genes.txt

    echo "Last but not least, let's merge with the desired phenotypes"
    # Replace spaces to prevent delimiter problems
    sed -i 's/ /_/g' ${PHENO}

    # Extract only WES-individuals
    head -1 ${PHENO} > ${TEMP}/WES_MCM_phenotypes.tab
    ${OVERLAP} ${FAM} 1 ${PHENO} 1 >> ${TEMP}/WES_MCM_phenotypes.tab
    
    # Merge with CMR-data
    ${MERGE} --file1 ${FULL_CMR} --file2 ${TEMP}/WES_MCM_phenotypes.tab --index f.eid > ${TEMP}/WES_MRI_MCM_phenotypes.txt

    ${MERGE} --file1 ${TEMP}/WES_MRI_MCM_phenotypes.txt --file2 ${TEMP}/${DIS}_IIDs_genes.txt --index f.eid | sed 's/ /\t/g' > ${OUT}/${DIS}_${SUF}

    echo "Have a great day!"
fi
