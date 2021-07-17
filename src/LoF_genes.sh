#!/bin/bash

head -1 /hpc/dhl_ec/aalasiri/scripts/WES_200K_LoF/data/loftk_WES_ukb_200K/WES_ukb_200K_v2/STAT_REPORT/Marion/LoF_gene_counts/WES_ukb_200K_ALLchrs_extr_gene.freq.filtered.counts | awk '{print "Column", $1, $2, $3, $4, "CM"}' > data/temp/LoF_genes.txt

while IFS= read -r line; do

  IFS=" " read -ra entries <<< "$line"
  GENE=$(echo "${entries[0]}")
  CM=$(echo "${entries[1]}")

  echo "Searching for the following term now: ${GENE}"
  echo -e "$(grep -wn ${GENE} /hpc/dhl_ec/aalasiri/scripts/WES_200K_LoF/data/loftk_WES_ukb_200K/WES_ukb_200K_v2/STAT_REPORT/Marion/LoF_gene_counts/WES_ukb_200K_ALLchrs_extr_gene.freq.filtered.counts | sed 's/:/\t/' | cut -f-5) ${CM}" >> data/temp/LoF_genes.txt
  echo ""

done < "data/temp/CM_genes.txt"
awk '$4 > 0 || $5 > 0 {print $3}' data/temp/LoF_genes.txt | sort -u > data/temp/LoF_genes_only.txt

cat /hpc/dhl_ec/aalasiri/scripts/WES_200K_LoF/data/loftk_WES_ukb_200K/WES_ukb_200K_v2/STAT_REPORT/Marion/LoF_gene_counts/WES_ukb_200K_ALLchrs_extr_gene.freq.filtered.counts | bin/transpose_perl.pl > data/temp/MCM_LoF_genes_transposed.txt

GNS=$(tail -n +2 data/temp/LoF_genes.txt | awk '$2 != "" {print $1}' | bin/transpose_perl.pl | sed 's/\t/,/g')

cut -f1,${GNS} data/temp/MCM_LoF_genes_transposed.txt | tail -n +2 > data/temp/MCM_LoF_genes.txt
sed -i 's/gene_symbol/f.eid/' data/temp/MCM_LoF_genes.txt

rm data/temp/LoF_genes.txt data/temp/LoF_genes_only.txt data/temp/MCM_LoF_genes_transposed.txt 