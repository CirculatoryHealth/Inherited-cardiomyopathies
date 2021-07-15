#!/bin/bash

plink=/hpc/local/CentOS7/dhl_ec/software/plink_v1.9


INPUT=/hpc/dhl_ec/aalasiri/lof/WES_ukb_50K # filtered 50 K WES
OUTPUT=/hpc/dhl_ec/aalasiri/Mendelian_CM_MRI_Exome/ARVC


sed 's/;/\t/g' WES_ARVC_mB.csv | cut -f3-4 | awk NF > snps_build38
sed 's/;/\t/g' WES_ARVC_mB.csv | cut -f1-4 | sed 's/ //g' | awk '$4 = " " {print $1,$2}' > snps_build37_has_no_build38
tail -n +2 snps_build37_has_no_build38 | awk -F"\t" '{print $1,$2}' | cut -d'-' -f1 |awk '$2 > 100 {print $0}'| sed 's/ /:/' | sort -u > snps_build37_has_no_build38.clean
tail -n +2 snps_build38 | awk -F"\t" '{print $1,$2}' | cut -d'-' -f1 | awk '$2 > 0 {print $0}' | sed 's/ /:/' | sort -u > snps_build38.clean # contains 429 SNPs

### convert built37 to 38
cut -d'-' -f1 remapped_remap_textarea | sort -u > converted_build37_to_build38 # contains 496 SNPs

## Merge both builds 
cat converted_build37_to_build38 > all_snps_build38.clean
cat snps_build38.clean >> all_snps_build38.clean
sort all_snps_build38.clean | uniq > all_snps_build38.clean1

## Get overlapped SNPs (build38 and converted build37)  between our list and all 200K UKB WES
/hpc/local/CentOS7/dhl_ec/software/overlap.pl all_snps_build38.clean.sorted 1 ukb_SNP_ID_WES_chrALL.txt_position 2 |  awk '{print $1}' > overlapped_ukb_ARVC_all_in_build38_snps # contains 227 SNPs
