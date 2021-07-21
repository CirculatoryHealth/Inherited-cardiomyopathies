#!/bin/bash
#
#SBATCH --job-name clean_pheno                                            		  # the name of this script
#SBATCH --output=/hpc/dhl_ec/mvanvugt/scripts/logs/clean_pheno.log              # log file of this script
#SBATCH --mem=10G                                                          		  # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#SBATCH --mail-user=m.vanvugt-2@umcutrecht.nl                                   # you can send yourself emails when the job is done; "--mail-user" and "--mail-type" go hand in hand
#SBATCH --mail-type=FAIL                                                        # you can choose: ALL=select all options; END=end of job; FAIL=abort of job


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
  echoerror "- Argument #1   --  Disease to perform analysis for, could be 'DCM'"
  echoerror "- Argument #2   --  Suffix of the output file, could be 'DCM'"
	echoerror ""
	echoerror "An example command would be: clean_pheno.sh [arg1: DCM]."
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  	# The wrong arguments are passed, so we'll exit the script now!
  	exit 1
}

if [[ $# -lt 2 ]]; then
  echo "Error, number of arguments found "$#"."
  script_arguments_error "You must supply [1] correct arguments when running this script"

else

  DIR="$1"
  ROOT=$(pwd)

  echo "Making phenotype file for the controls"
  bin/merge_tables.pl --file1 data/temp/WES_MRI_MCM_phenotypes.txt --file2 ${DIR}/Control_IDs.txt --index f.eid > data/raw/Controls_full.txt
  echo ""
  echo "Combining the phenotype files"
  Rscript --vanilla ${ROOT}/src/Combine_pheno.R data/temp/ _raw.txt MCM_raw_full
  echo ""
  echo "Summarizing genetic information"
  Rscript --vanilla ${ROOT}/src/MCM_gene_summary.R data/raw/MCM_clean_full.rds data/raw data/processed/ MCM
  echo ""
  echo "Cleaning up the phenotype file"
  Rscript --vanilla ${ROOT}/src/MCM_pheno_clean.R data/raw/MCM_clean_full.rds data/raw data/processed MCM
  echo ""
  echo "Merging the cleaned files"
  bin/merge_tables.pl --file1 data/processed/MCM_gene_summary.tsv --file2 data/processed/MCM_cleaned.tsv --index f.eid | sed 's/ /\t/g' > data/processed/MCM_final_pheno.tsv

fi
