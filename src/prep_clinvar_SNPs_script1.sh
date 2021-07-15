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
	echoerror ""
	echoerror "An example command would be: prep_clinvar_SNPs_script1.sh [arg1]"
	echoerror "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
 	echo ""
	script_copyright_message
	exit 1
}

echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "+            Extraction of [Likely] pathogenic mutations for cardiomyopathies from Clinvar              +"
echobold "+                                                                                                       +"
echobold "+                                                                                                       +"
echobold "+ * Written by  : Abdulrahman Alasiri                                                                   +"
echobold "+ * E-mail      : a.i.alasiri@umcutrecht.nl                                                             +"
echobold "+ * Last update : 2021-02-22                                                                            +"
echobold "+ * Version     : 0.1.0                                                                                 +"
echobold "+                                                                                                       +"
echobold "+ * Description : This script will extract [Likely] pathogenic mutations                                +"
echobold "+                 for ARVC, DCM and HCM from Clinvar.                                                   +"
echobold "+                                                                                                       +"
echobold "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's date and time: "$(date)
TODAY=$(date +"%Y%m%d")
echo ""

### START of if-else statement for the number of command-line arguments passed ###
if [[ $# -lt 1 ]]; then
    echoerrorflash "                                     *** Oh no! Computer says no! ***"
    echo ""
    script_arguments_error "You must supply at least [1] argument when running extraction of cardiomyopathy mutations from Clinvar!"

else
    ### REQUIRED | GENERALS
    DIS="$1"    # Disease (DCM/HCM/ARVC)

    ### TOOLS
    OVERLAP=/hpc/local/CentOS7/dhl_ec/software/overlap.pl
    MERGE=/hpc/local/CentOS7/dhl_ec/software/merge_tables.pl

    ### PROJECT SPECIFIC
    PROJECTDIR=/hpc/dhl_ec/aalasiri/scripts/Inherited_cardiomyopathies # the root directoryto repo
    INPUT=${PROJECTDIR}/data/raw
    OUTPUT=${PROJECTDIR}/data/processed
    TEMP=${PROJECTDIR}/data/temp
    UKB_200K_WES=/hpc/ukbiobank/WES_200K_2020


### Directories
#MAIN_DIR=/hpc/dhl_ec/aalasiri/Mendelian_CM_MRI_Exome_200K
#DIR=${MAIN_DIR}/${DIS}
#UKB_200K_WES=/hpc/ukbiobank/WES_200K_2020


### STEP 1 ###
## Extract snps loci  from ClinVar
    echo "ID GRCh38Position" > ${TEMP}/${DIS}_P_and_LP_positionID
    cat ${INPUT}/${DIS}_P_and_LP_v3.txt | sed 's/ /_/g6' | awk -F"\t" '{print $15,$10,$11,$15}' | cut -d ' ' -f1,2,3 | sed 's/:/ /1; s/ /:/3; s/_-_/ /g' | awk '{print $2,$3,$2}' | sed 's/:/ /4' | awk '{print $1,$2":"$4}' | tail -n +2 >> ${TEMP}/${DIS}_P_and_LP_positionID

## SNP position ID for UKB WES 200K (NO NEED to do it if the file already created)
    if [ ! -e ${INPUT}/ukb_SNP_ID_WES_chrALL.txt ]; then
      for CHR in $(seq 1 22); do
        awk '{print $2, $1":"$4":"$6":"$5}' ${UKB_200K_WES}/UKBexomeOQFE_200K_chr${CHR}_v1.bim >> ${TEMP}/${DIS}/ukb_SNP_ID_WES_chrALL.txt
      done
    fi

## Overlap with WES
    ${OVERLAP} ${TEMP}/${DIS}_P_and_LP_positionID 2 ${INPUT}/ukb_SNP_ID_WES_chrALL.txt 2 | cut -d' ' -f1 >  ${OUTPUT}/${DIS}_overlap_path_SNPs_WES_SNPs.txt
    sed 's/:/ /g5' ${INPUT}/ukb_SNP_ID_WES_chrALL.txt > ${TEMP}/${DIS}_ukb_SNP_ID_WES_chrALL.txt_position # replace ":" between position and ref to space
    sed 's/:/ /g4' ${DIR}/${DIS}_P_and_LP_positionID > ${TEMP}/${DIS}_P_and_LP_positionID_position
    ${OVERLAP} ${TEMP}/${DIS}_P_and_LP_positionID_position 2 ${TEMP}/${DIS}_ukb_SNP_ID_WES_chrALL.txt_position 2 > ${OUTPUT}/${DIS}_overlap_path_SNPs_WES_SNPs.txt_position # get overlap by position

### In screen 1
##/hpc/local/CentOS7/dhl_ec/software/overlap.pl ${OUTPUT}/${DIS}_overlap_path_SNPs_WES_SNPs.txt 1 ${OUTPUT}/${DIS}_overlap_path_SNPs_WES_SNPs.txt_position 1 -v | awk '$1 ~ /I/ {print $0}' | less # get non-overlap Insersion
##/hpc/local/CentOS7/dhl_ec/software/overlap.pl ${OUTPUT}/${DIS}_overlap_path_SNPs_WES_SNPs.txt 1 ${OUTPUT}/${DIS}_overlap_path_SNPs_WES_SNPs.txt_position 1 -v | awk '$1 ~ /D/ {print $0}' | less  # get non-overlap Deletion
### in screen 2
##less ${DIS}_P_and_LP_positionID ## search in 'less' usin "/" then paste each position from screen 1 "column2", then check if the indel is the same
## compare indel in overlap to the one in ${DIS}_P_and_LP_positionID file, then write them in excel sheet
## emacs overlap_path_SNPs_WES_SNPs.txt # add matched indel snps

### FOR ARVC, combine all SNPs in these two files:
## cat overlapped_ukb_ARVC_all_in_build38_snps >> overlap_path_SNPs_WES_SNPs.txt

### After completing this part, run "run_clinvarSNPs_WES_UKB_script2.sh"
