#!/bin/bash
### PROJECT SPECIFIC
PROJECTDIR=/hpc/dhl_ec/aalasiri/scripts/Inherited_cardiomyopathies # the root directoryto repo
INPUT=${PROJECTDIR}/data/raw
OUTPUT=${PROJECTDIR}/data/processed
TEMP=${PROJECTDIR}/data/temp
BIN=${PROJECTDIR}/bin


  sed 's/;/\t/g' ${INPUT}/WES_ARVC_mB.csv | cut -f3-4 | awk NF > ${TEMP}/snps_build38
sed 's/;/\t/g' ${INPUT}/WES_ARVC_mB.csv | cut -f1-4 | sed 's/ //g' | awk '$4 = " " {print $1,$2}' > ${TEMP}/snps_build37_has_no_build38
tail -n +2 ${TEMP}/snps_build37_has_no_build38 | awk -F"\t" '{print $1,$2}' | cut -d'-' -f1 |awk '$2 > 100 {print $0}'| sed 's/ /:/' | sort -u > ${TEMP}/snps_build37_has_no_build38.clean
tail -n +2 ${TEMP}/snps_build38 | awk -F"\t" '{print $1,$2}' | cut -d'-' -f1 | awk '$2 > 0 {print $0}' | sed 's/ /:/' | sort -u > ${TEMP}/snps_build38.clean # contains 429 SNPs

### convert built37 to 38
cut -d'-' -f1 ${BIN}/remapped_remap_textarea | sort -u > ${TEMP}/converted_build37_to_build38 # contains 496 SNPs

## Merge both builds
cat ${TEMP}/converted_build37_to_build38 > ${TEMP}/ARVC_all_snps_build38.clean
cat ${TEMP}/snps_build38.clean >> ${TEMP}/ARVC_all_snps_build38.clean
sort -n -k1.1 ${TEMP}/ARVC_all_snps_build38.clean | uniq > ${TEMP}/ARVC_all_snps_build38.clean.sorted


## Get overlapped SNPs (build38 and converted build37)  between our list and all 200K UKB WES
/hpc/local/CentOS7/dhl_ec/software/overlap.pl ${TEMP}/ARVC_all_snps_build38.clean.sorted 1 ${TEMP}/ARVC_ukb_SNP_ID_WES_chrALL.txt_position |  awk '{print $1}' > ${OUTPUT}/overlapped_ukb_ARVC_all_in_build38_snps # contains 227 SNPs
## combine all SNPs in these two files:
cat ${OUTPUT}/overlapped_ukb_ARVC_all_in_build38_snps >> ${OUTPUT}/ARVC_overlap_path_SNPs_WES_SNPs.txt
