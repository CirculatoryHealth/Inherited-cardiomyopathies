#!/bin/bash

ROOT=/hpc/dhl_ec/mvanvugt/UKBB/Inherited-cardiomyopathies

echo "Getting overlapping individuals"
${ROOT}/bin/overlap.pl ${ROOT}/data/temp/DCM_temp/DCM_gene.vcf 1 ${ROOT}/data/temp/HCM_temp/HCM_gene.vcf 1 | cut -f1,2 | tail -n +2 > ${ROOT}/data/temp/overlap_DCM_HCM.txt
echo "f.eid SNP Gene" > ${ROOT}/data/temp/ID_overlap_DCM_HCM.txt
while IFS= read -r line; do

    entries=($line)
    SNP=$(echo "${entries[0]}")
    GENE=$(echo "${entries[1]}")

    echo "${SNP}"
    NUM=`grep ${SNP} ${ROOT}/data/temp/HCM_temp/HCM_gene.vcf | ${ROOT}/bin/transpose_perl.pl | grep -n 1 | tail -n +2 | sed 's/:/\t/g' | cut -f1 | ${ROOT}/bin/transpose_perl.pl | sed 's/\t/,/g'`
    head -1 ${ROOT}/data/temp/HCM_temp/HCM_gene.vcf | cut -f${NUM} | ${ROOT}/bin/transpose_perl.pl | awk -v snp="${SNP}" -v gen="${GENE}" '{print $1, snp, gen}' >> ${ROOT}/data/temp/ID_overlap_DCM_HCM.txt

done < "${ROOT}/data/temp/overlap_DCM_HCM.txt"

${ROOT}/bin/overlap.pl ${ROOT}/data/temp/DCM_temp/DCM_gene.vcf 1 ${ROOT}/data/temp/ACM_temp/ACM_gene.vcf 1 | cut -f1,2 | tail -n +2 > ${ROOT}/data/temp/overlap_DCM_ACM.txt
echo "f.eid SNP Gene" > ${ROOT}/data/temp/ID_overlap_DCM_ACM.txt
while IFS= read -r line; do

    entries=($line)
    SNP=$(echo "${entries[0]}")
    GENE=$(echo "${entries[1]}")

    echo "${SNP}"
    NUM=`grep ${SNP} ${ROOT}/data/temp/DCM_temp/DCM_gene.vcf | ${ROOT}/bin/transpose_perl.pl | grep -n 1 | tail -n +2 | sed 's/:/\t/g' | cut -f1 | ${ROOT}/bin/transpose_perl.pl | sed 's/\t/,/g'`
    head -1 ${ROOT}/data/temp/DCM_temp/DCM_gene.vcf | cut -f${NUM} | ${ROOT}/bin/transpose_perl.pl | awk -v snp="${SNP}" -v gen="${GENE}" '{print $1, snp, gen}' >> ${ROOT}/data/temp/ID_overlap_DCM_ACM.txt

done < "${ROOT}/data/temp/overlap_DCM_ACM.txt"

echo "Merge with phenotypes"
${ROOT}/bin/merge_tables.pl --file1 ${ROOT}/data/processed/MCM_final_pheno.tsv --file2 ${ROOT}/data/temp/ID_overlap_DCM_ACM.txt --index f.eid | sed 's/ /\t/g' > ${ROOT}/data/processed/overlap_DCM_ACM_pheno.tsv

${ROOT}/bin/merge_tables.pl --file1 ${ROOT}/data/processed/MCM_final_pheno.tsv --file2 ${ROOT}/data/temp/ID_overlap_DCM_HCM.txt --index f.eid | sed 's/ /\t/g' > ${ROOT}/data/processed/overlap_DCM_HCM_pheno.tsv
